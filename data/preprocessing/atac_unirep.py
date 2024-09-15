r"""
Compute unimodal representation of ATAC data
(for computing the neighbor conservation metric)
"""

import argparse
import pathlib

import anndata
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
scmbench_dir = os.path.abspath(os.path.join(parent_dir,'..'))
print('scmbench_dir',scmbench_dir)
sys.path.append(scmbench_dir)
import SCMBench


def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Compute unimodal representation of ATAC data"
    )
    parser.add_argument(
        "-i", "--input", dest="input", type=pathlib.Path, required=True,
        help="Path to input file (.h5ad)"
    )
    parser.add_argument(
        "-d", "--dim", dest="dim", type=int, default=50,
        help="Dimensionality of the representation"
    )
    parser.add_argument(
        "-o", "--output", dest="output", type=pathlib.Path, required=True,
        help="Path to output file (.h5ad)"
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    print("[1/3] Reading data...")
    adata = anndata.read_h5ad(args.input)
    print("[2/3] Computing LSI...")
    SCMBench.data.lsi(adata, n_components=args.dim, use_highly_variable=False, n_iter=15)
    print("[3/3] Saving results...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    anndata.AnnData(
        X=adata.obsm["X_lsi"], obs=adata.obs
    ).write(args.output, compression="gzip")


if __name__ == "__main__":
    main(parse_args())
