#! /usr/bin/env python3
# Author: Kristopher Kieft
# University of Wisconsin-Madison
# 2022

# modified by cody to fix an issue


import os
import numpy as np

np.seterr(all="raise")
import argparse
from fasta_parse import fasta_parse
from numba import jit


def arguments():
    descript = """
        Calculate average coverage per scaffold from samtools depth output (or similar format)
        with a single or multiple samples.

        Coverage file: (include header, or --no_header)
        Col 1: scaffold
        Col 2: position (1-based)
        Col 3-n: coverage per sample

        Regions of interest file:
        Col 1: scaffold
        Col 2: name of region
        Col 3: start coordinate (1-based)
        Col 4: end coordinate (1-based)

        Usage:
        cov_by_region.py -i coverage.tsv -r regions.tsv -o output.tsv -f scaffolds.fasta
    """
    parser = argparse.ArgumentParser(
        description=descript,
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-i",
        metavar="",
        type=str,
        nargs=1,
        required=True,
        help="input coverage file (include header)",
    )
    parser.add_argument(
        "-f",
        metavar="",
        type=str,
        nargs=1,
        required=True,
        help="input fasta file of represented sequences (for seq length)",
    )
    parser.add_argument(
        "-r",
        metavar="",
        type=str,
        nargs=1,
        required=True,
        help="regions of interest tsv file",
    )
    parser.add_argument(
        "-c",
        metavar="",
        type=str,
        help="Vibrant prophage coordinates file (it's ok if it has multiple headers). This is optional if you manually shift your prophage coordinates to start from 1 before using this script.",
    )
    parser.add_argument(
        "-o",
        metavar="",
        type=str,
        nargs=1,
        required=True,
        help="output coverage information tsv",
    )
    parser.add_argument(
        "--no_header",
        action="store_true",
        help="coverage file has no header and is single sample",
    )

    args = parser.parse_args()
    covfile = args.i[0]
    fasta = args.f[0]
    regionfile = args.r[0]
    coordsfile = args.c
    outfile = args.o[0]

    if os.path.exists(outfile):
        print("\nError: output file already exists. Exiting.\n")
        exit()

    return covfile, regionfile, coordsfile, outfile, fasta, args.no_header


def get_samples(f):
    with open(f, "r") as infile:
        header = infile.readline().strip("\n").split("\t")
        samples = len(header) - 2
    return samples


def get_lengths(fasta):
    lengths = {}
    for name, seq in fasta_parse(fasta):
        seq = seq.replace("\n", "")
        lengths[name] = len(seq)

    return lengths


def get_regions(regionfile, coords):
    """regions file has these columns:

      - scaffold
      - region, basically a unique ID for a subregion in the scaffold
      - start
      - stop
    """
    regions = {}
    with open(regionfile, "r") as infile:
        next(infile)  # skip header im guessing
        for line in infile:
            line = line.rstrip().split("\t")
            # genome: (region, start, end) <- idk what this means

            # this is just the scaffold name
            genome = line[0]

            # region ID
            r = line[1]
            start = int(line[2])
            end = int(line[3])

            if "_fragment_" in genome:
                # edit start/stop of "fragments" only
                # fragments are prophages
                # vibrant coordinates are 1 based
                partial_name = genome.split("__")[-1]
                fragment_start = coords.get(partial_name, 0)
                # if it is a prophage shift the start
                # to be relatively at 1
                start -= fragment_start
                end -= fragment_start
            else:
                # shift everything else to 0 based coords?
                # shouldn't we do this with prophages too??
                # i believe this is an error to be in the else block
                # so I am removing the else block since I think all
                # coords need to be shifted by 1 back to get 0-indexed
                # coordinates for the coverage file
                # according to the description, the coverage file
                # is also 1 based??
                start -= 1

            regions.setdefault(genome, []).append((r, start, end))  # now zero based

    return regions


def get_coords(coordsfile):
    """coords file from VIBRANT has these columns

      - 0 scaffold
      - 1 fragment
      - 2 protein start
      - 3 protein stop
      - 4 protein length
      - 5 nucleotide start
      - 6 nucleotide stop
      - 7 nucleotide length

    Thus `line[1]` is the fragment name, ie the ID for an excised fragment from a scaffold. `line[5]` is the nucleotide start.

    Returns:
        dict[str, int]: maps fragment name to nucleotide stop coordinate
    """
    coords = {}
    with open(coordsfile, "r") as f:
        for line in f:
            line = line.strip("\n").split("\t")
            # fragment: (start, end)
            if line[1] != "fragment":  # header
                coords[line[1]] = int(line[5])

    return coords


@jit(nopython=True)
def add_depth(depth, pos, coverages):
    depth[:, pos - 1] = coverages
    return depth


def get_avg(depth, regions, hold_scaff, samples):
    avg_scaffold, avg_scaffold_partial, avg_regions = (
        [],
        [],
        {},
    )  # np arr, np arr, list of (region, np arr)

    depth_partial = depth.copy()
    avg_scaffold = np.mean(depth, axis=1)  # whole scaffold

    for entry in regions.get(hold_scaff, []):
        region, start, stop = entry
        for i in range(samples):
            avg_regions.setdefault(region, []).append(np.mean(depth[i][start:stop]))
            for n in range(start, stop):
                depth_partial[i][n] = np.nan

    avg_scaffold_partial = np.nanmean(
        depth_partial, axis=1
    )  # whole scaffold excluding regions

    return avg_scaffold, avg_scaffold_partial, avg_regions


def make_header(samples):
    w = "\t".join([f"avg whole scaffold #{i+1}" for i in range(samples)])
    p = "\t".join([f"avg partial scaffold #{i+1}" for i in range(samples)])
    r = "\t".join([f"avg region #{i+1}" for i in range(samples)])
    header = f"scaffold\tregion\t{w}\t{p}\t{r}\n"

    return header


def main(covfile, regionfile, coordsfile, outfile, fasta, noheader):
    if not noheader:
        samples = get_samples(covfile)
    else:
        samples = 1
    header = make_header(samples)
    lengths = get_lengths(fasta)

    # if coords file is None just use an empty dict
    coords = get_coords(coordsfile) if coordsfile is not None else dict()

    # if coords file is None, the default shift will be 0,
    # since if the coords file is not provided, that means the
    # the prophage coordinates were already adjusted
    regions = get_regions(regionfile, coords)

    with open(covfile, "r") as infile, open(outfile, "w") as out:
        out.write(header)
        for line in infile:
            line = line.strip("\n").split("\t")
            scaffold = line[0]
            position = int(line[1])
            coverages = np.array([c for c in line[2:]], int)

            try:
                if scaffold == hold_scaff:
                    depth = add_depth(depth, position, coverages)

                else:
                    avg_scaffold, avg_scaffold_partial, avg_regions = get_avg(
                        depth, regions, hold_scaff, samples
                    )
                    avg_scaffold = "\t".join([str(a) for a in avg_scaffold])
                    avg_scaffold_partial = "\t".join(
                        [str(a) for a in avg_scaffold_partial]
                    )

                    if avg_regions:
                        for entry, avg_values in avg_regions.items():
                            # r, avg_entry = entry
                            avg_values = "\t".join([str(a) for a in avg_values])
                            out.write(
                                f"{hold_scaff}\t{entry}\t{avg_scaffold}\t{avg_scaffold_partial}\t{avg_values}\n"
                            )
                    else:
                        out.write(f"{hold_scaff}\tNA\t{avg_scaffold}\tNA\tNA\n")

                    del lengths[hold_scaff]
                    length = lengths[scaffold]
                    depth = np.zeros((samples, length))
                    depth = add_depth(depth, position, coverages)

            except NameError:
                length = lengths[scaffold]
                depth = np.zeros((samples, length))
                depth = add_depth(depth, position, coverages)

            hold_scaff = scaffold

        # last scaffold
        avg_scaffold, avg_scaffold_partial, avg_regions = get_avg(
            depth, regions, hold_scaff, samples
        )
        avg_scaffold = "\t".join([str(a) for a in avg_scaffold])
        avg_scaffold_partial = "\t".join([str(a) for a in avg_scaffold_partial])

        if avg_regions:
            for entry, avg_values in avg_regions.items():
                # r, avg_entry = entry
                avg_values = "\t".join([str(a) for a in avg_values])
                out.write(
                    f"{hold_scaff}\t{entry}\t{avg_scaffold}\t{avg_scaffold_partial}\t{avg_values}\n"
                )
        else:
            out.write(f"{hold_scaff}\tNA\t{avg_scaffold}\tNA\tNA\n")

        # no coverage
        for scaff in lengths.keys():
            out.write(f"{scaff}\tNA\t0\tNA\tNA\n")


if __name__ == "__main__":
    covfile, regionfile, coordsfile, outfile, fasta, noheader = arguments()
    main(covfile, regionfile, coordsfile, outfile, fasta, noheader)
