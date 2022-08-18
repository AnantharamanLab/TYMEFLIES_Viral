#!/usr/bin/env python3
import argparse
import multiprocessing
from itertools import repeat
from typing import Iterator, List, Tuple

import numpy as np
from pandas import to_numeric
from sklearn.linear_model import LinearRegression

_HEADER = [
    "var",
    "reg_slope",
    "reg_yint",
    "reg_rsquared",
]

NULLDATA = "\t".join("NA" for _ in range(len(_HEADER) - 2))


def read_data(file: str) -> Tuple[str, np.ndarray]:
    """Read tab-delimited file that has 1 row as a single variable that will
    be regressed against the index position.

    Args:
        file (str): data file as described above

    Returns:
        Tuple[str, np.ndarray]: (var_name, data_array)
    """
    with open(file) as fp:
        var, *values = fp.readline().split("\t")

    array = to_numeric(np.array(values), errors="coerce")
    return var, array[~np.isnan(array)]


def regression(array: np.ndarray) -> Tuple[float, float, float]:
    """Perform a simple linear regression for the variable in `array`
    and return the slope, y-intercept, and R-squared value regressed
    against the index position.

    Args:
        array (np.ndarray): 1 dimension array

    Returns:
        Tuple[float, float, float]: (slope, y-int, R-sq)
    """
    X: np.ndarray = array.reshape(-1, 1)
    y: np.ndarray = np.arange(start=1, stop=array.size + 1).reshape(-1, 1)
    model: LinearRegression = LinearRegression().fit(X, y)
    rsquare: float = model.score(X, y)
    slope: float = model.coef_[0][0]
    yint: float = model.intercept_[0]
    return slope, yint, rsquare


def main(file: str, output: str, min_number: int) -> None:
    HEADER = "\t".join(_HEADER)
    with open(output, "w") as outfile:
        outfile.write(f"{HEADER}\n")
        var, array = read_data(file)
        nelems = array.size
        if nelems < min_number:
            outfile.write(f"{var}\t{NULLDATA}\n")
        else:
            slope, yint, rsquare = regression(array)
            outfile.write(f"{var}\t{slope:.4f}\t{yint:.4f}\t{rsquare:.4f}\n")


def read_batch_file(file: str) -> Iterator[Tuple[str, str]]:
    with open(file) as fp:
        for line in fp:
            infile, outfile = line.rstrip().split("\t")
            yield infile, outfile


def _parallel_helper(files: List[str], min_number: int) -> None:
    file, output = files
    main(file, output, min_number)


def main_parallel(file: str, min_number: int, jobs: int) -> None:
    files = read_batch_file(file)
    with multiprocessing.Pool(processes=jobs) as pool:
        pool.starmap(_parallel_helper, zip(files, repeat(min_number)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the spearman correlation between two rows of data",
    )

    output_args = parser.add_mutually_exclusive_group(required=True)

    input_help = (
        "SINGLE FILE: tab-delimited file with a single row of data where the first column is the variable name. 'NA' is a null value. Requires -o argument. "
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
        help="name of output tab-delimited file that will contain the two variable names, the correlation coefficient, and the pvalue",
    )
    output_args.add_argument(
        "--batch-mode",
        default=False,
        action="store_true",
        help="use if providing a batch tab-delimited file with the input files on the left and output paths on the right",
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
        main_parallel(args.input, args.num_limit, args.jobs)
    else:
        main(args.input, args.output, args.num_limit)
