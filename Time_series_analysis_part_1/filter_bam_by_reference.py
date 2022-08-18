#!/usr/bin/env python3
DESCRIPTION = """
Author:       Cody Martin

Affilitation: University of Wisconsin-Madison
              Department of Bacteriology
              Anantharaman lab

Purpose: Given a text file of scaffold names to keep,
         filter a bam file to keep only the alignments
         that map to a given reference in the scaffolds.
         Additionally, only keep the sequence headers
         corresponding to your scaffolds of interest.

Usage:

1) Filter one bam file
python3 filter_bam_by_reference.py -b <in.bam> -r <scaffolds.txt>

2a) Filter multiple bam files
python3 filter_bam_by_reference.py -b <in1.bam> <in2.bam> -r <scaffolds.txt>

2b) Filter multiple bam files with unix wildcard pattern
python3 filter_bam_by_reference.py -b ./*.bam -r <scaffolds.txt>
"""
import pysam
import argparse
import multiprocessing
from functools import partial
from typing import List


def read_references(file: str) -> set:
    """Read reference names into a set for constant time
    membership lookup.

    Args:
        file (str): file name of the newline-delimited scaffolds

    Returns:
        set: converts text file of reference names into a set
    """
    with open(file) as fp:
        references = {line.rstrip() for line in fp}
    return references


def filter_header(
    bamfile: pysam.AlignmentFile, scaffolds: set
) -> pysam.AlignmentHeader:
    """Filter the header of the BAM file by the set of `scaffolds` to keep.

    Args:
        bamfile (pysam.AlignmentFile): name of BAM file
        scaffolds (set): set of reference scaffold names to keep

    Returns:
        pysam.AlignmentHeader: filtered header object
    """
    #### Filter the SAM header ####
    #### otherwise, the filtered BAM files would all still
    #### have the original large header.
    filtered_header = dict()
    for tag, metadata in bamfile.header.items():
        if tag == "SQ":
            filtered_header["SQ"] = [
                reference for reference in metadata if reference["SN"] in scaffolds
            ]
        else:
            filtered_header[tag] = metadata

    return pysam.AlignmentHeader.from_dict(filtered_header)


def filter_bam(bamfile: str, scaffolds: set) -> None:
    """Filter a BAM file with reference names present in the
    `scaffolds` set. This will also filter the header.

    Args:
        bamfile (str): name of a BAM file
        scaffolds (set): set of reference scaffold names to keep
    """
    with pysam.AlignmentFile(bamfile, "rb") as inbam:
        filtered_header = filter_header(inbam, scaffolds)

        #### Filter reads ####
        output = f'{bamfile.split(".bam")[0]}.filtered.bam'
        with pysam.AlignmentFile(output, "wb", header=filtered_header) as outbam:
            for alignment in inbam.fetch(until_eof=True):
                if alignment.reference_name in scaffolds:
                    # need reference id numbers to match up
                    # with the new filtered header
                    # reference_id and next_reference_id both are updated to keep paired read info similar
                    alignment.reference_id = filtered_header.get_tid(
                        alignment.reference_name
                    )
                    alignment.next_reference_id = filtered_header.get_tid(
                        alignment.reference_name
                    )
                    outbam.write(alignment)


def main(bamfiles: List[str], referencefile: str, jobs: int) -> None:
    """Filter BAM files in parallel.

    Args:
        bamfiles (List[str]): list of BAM file names
        referencefile (str): text file of reference / scaffold names to keep
        jobs (int): number of parallel jobs to run
    """
    scaffolds = read_references(referencefile)

    # don't make more processes than needed
    if jobs > len(bamfiles):
        jobs = len(bamfiles)

    filter_bam_by_references = partial(filter_bam, scaffolds=scaffolds)

    with multiprocessing.Pool(processes=jobs) as pool:
        pool.map(filter_bam_by_references, bamfiles)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-b",
        "--bamfiles",
        metavar="BAM",
        nargs="+",
        required=True,
        help="bam files to filter. This can be a unix-style wildcard (*.bam) or a single file name.",
    )

    parser.add_argument(
        "-r",
        "--references",
        required=True,
        help="text file of reference names to keep. Each reference must be on a newline",
    )

    parser.add_argument(
        "-j",
        "--jobs",
        default=5,
        type=int,
        help="max number of parallel jobs to run, (default: %(default)s)",
    )

    args = parser.parse_args()
    bamfiles = args.bamfiles
    referencefile = args.references
    jobs = args.jobs
    main(bamfiles, referencefile, jobs)
