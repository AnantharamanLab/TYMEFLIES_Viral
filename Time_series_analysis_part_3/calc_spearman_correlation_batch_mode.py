#!/usr/bin/env python3
from itertools import repeat
from typing import Iterator, List, Tuple
import pandas as pd
from scipy.stats import spearmanr
import argparse
import multiprocessing

_HEADER = [
    "var1",
    "var2",
    "spearman_pval",
    "spearman_corr",
]

NULLDATA = "\t".join("NA" for _ in range(len(_HEADER) - 2))


def read_data(file: str, remove_zeros: bool) -> pd.DataFrame:
    """Read tab-delimited file that has exactly 2 rows corresponding to
    2 different variables. The columns are values from each sample that will
    be regressed and correlated to each other.

    Args:
        file (str): data file as described above
        remove_zeros (bool): remove columns with any 0.0 values

    Returns:
        pd.DataFrame: data read into a pandas dataframe
    """
    data = (
        pd.read_table(file, header=None)
        .rename({0: "Variable"}, axis=1)
        .set_index("Variable")
        .dropna(axis=1)
    )
    if remove_zeros:
        return data.drop(columns=data.columns[(data == 0.0).any()])
    return data


def spearman_correlation(df: pd.DataFrame) -> Tuple[float, float]:
    """Calculate the spearmean correlation using `scipy.stats.spearmanr`.

    Args:
        df (pd.DataFrame): 2xN dataframe, where the spearman correlation
            will be calculated on the Nx2 transposed dataframe

    Returns:
        Tuple[float, float]: (pvalue, correlation)
    """
    correlation, pvalue = spearmanr(df.T)
    return pvalue, correlation


def main(file: str, output: str, min_number: int, remove_zeros: bool) -> None:
    HEADER = "\t".join(_HEADER)
    with open(output, "w") as outfile:
        outfile.write(f"{HEADER}\n")
        data = read_data(file, remove_zeros)
        ncols = data.shape[1]
        variables = "\t".join(data.index.to_list())
        if ncols < min_number:
            outfile.write(f"{variables}\t{NULLDATA}\n")
        else:
            pvalue, corr = spearman_correlation(data)
            outfile.write(f"{variables}\t{pvalue:.4f}\t{corr:.4f}\n")


def read_batch_file(file: str) -> Iterator[Tuple[str, str]]:
    with open(file) as fp:
        for line in fp:
            infile, outfile = line.rstrip().split("\t")
            yield infile, outfile


def _parallel_helper(files: List[str], min_number: int, remove_zeros: bool) -> None:
    file, output = files
    main(file, output, min_number, remove_zeros)


def main_parallel(file: str, min_number: int, remove_zeros: bool, jobs: int):
    files = read_batch_file(file)
    with multiprocessing.Pool(processes=jobs) as pool:
        pool.starmap(
            _parallel_helper, zip(files, repeat(min_number), repeat(remove_zeros))
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the spearman correlation between two rows of data"
    )

    output_args = parser.add_mutually_exclusive_group(required=True)

    input_help = (
        "SINGLE FILE: tab-delimited file with two row of data where the first column is the variable name. 'NA' is a null value. Requires -o argument. "
        "BATCH MODE: input file paths and output file paths for each individual file to be processed that are in the valid format listed above. These two paths are separated by a tab. Requires --batch-mode flag."
    )

    parser.add_argument(
        "-i", "--input", required=True, help=input_help,
    )
    parser.add_argument(
        "-n",
        "--num-limit",
        default=5,
        help="minimum number of columns remaining to do processing (default: %(default)s)",
    )
    output_args.add_argument(
        "-o",
        "--output",
        help="name of output tab-delimited file that will contain the two variable names, the pvalue, and the correlation coefficient",
    )
    output_args.add_argument(
        "--batch-mode",
        default=False,
        action="store_true",
        help="use if providing a batch tab-delimited file with the input files on the left and output paths on the right",
    )
    parser.add_argument(
        "--remove-zeros",
        default=False,
        action="store_true",
        help="use to remove columns that contain 0 (default: %(default)s)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=20,
        help="number of parallel jobs to run (default: %(default)s)",
    )
    args = parser.parse_args()
    if args.batch_mode:
        main_parallel(args.input, args.num_limit, args.remove_zeros, args.jobs)
    else:
        main(args.input, args.output, args.num_limit, args.remove_zeros)
