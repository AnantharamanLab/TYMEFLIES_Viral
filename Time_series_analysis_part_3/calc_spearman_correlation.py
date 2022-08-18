#!/usr/bin/env python3
from typing import Tuple
import pandas as pd
from scipy.stats import spearmanr
import argparse

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
        print(data)
        ncols = data.shape[1]
        variables = "\t".join(data.index.to_list())
        if ncols < min_number:
            outfile.write(f"{variables}\t{NULLDATA}\n")
        else:
            pvalue, corr = spearman_correlation(data)
            outfile.write(f"{variables}\t{pvalue:.4f}\t{corr:.4f}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the spearman correlation between two rows of data"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="input tab-delimited file that contains only two rows representing the variables, 'NA' as the only null values",
    )
    parser.add_argument(
        "-n",
        "--num-limit",
        default=5,
        help="minimum number of columns remaining to do processing (default: %(default)s)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="name of output tab-delimited file that will contain the two variable names, the pvalue, and the correlation coefficient",
    )
    parser.add_argument(
        "--remove-zeros",
        default=False,
        action="store_true",
        help="use to remove columns that contain 0 (default: %(default)s)",
    )
    args = parser.parse_args()
    main(args.input, args.output, args.num_limit, args.remove_zeros)
