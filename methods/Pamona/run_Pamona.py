#!/usr/bin/env python

r"""
Run UnionCom
"""

import argparse
import pathlib
import sys
import time

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

import SCMBench

sys.path.append("../custom")
import pamona
from pamona import Pamona


def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-rna", dest="input_rna", type=pathlib.Path, required=True,
        help="Path to input RNA dataset (.h5ad)"
    )
    parser.add_argument(
        "--input-atac", dest="input_atac", type=pathlib.Path, required=True,
        help="Path to input ATAC dataset (.h5ad)"
    )
    parser.add_argument(
        "-s", "--random-seed", dest="random_seed", type=int, default=0,
        help="Random seed"
    )
    parser.add_argument(
        "--output-rna", dest="output_rna", type=pathlib.Path, required=True,
        help="Path of output RNA latent file (.csv)"
    )
    parser.add_argument(
        "--output-atac", dest="output_atac", type=pathlib.Path, required=True,
        help="Path of output ATAC latent file (.csv)"
    )
    parser.add_argument(
        "-r", "--run-info", dest="run_info", type=pathlib.Path, required=True,
        help="Path of output run info file (.yaml)"
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    np.random.seed(args.random_seed)

    print("[1/4] Reading data...")
    rna = anndata.read_h5ad(args.input_rna)
    atac = anndata.read_h5ad(args.input_atac)

    print("[2/4] Preprocessing...")
    start_time = time.time()
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna, max_value=10)
    sc.tl.pca(
        rna, n_comps=min(100, rna.shape[0]),
        use_highly_variable=True, svd_solver="auto"
    )
    SCMBench.data.lsi(
        atac, n_components=min(100, atac.shape[0]),
        use_highly_variable=False, n_iter=15
    )

    print("[3/4] Training Pamona...")
    output_dir = args.run_info.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    X = rna.obsm["X_pca"]
    Y = atac.obsm["X_lsi"]
    X = pamona.utils.zscore_standardize(np.asarray(X))
    Y = pamona.utils.zscore_standardize(np.asarray(Y))
    Pa = Pamona.Pamona(manual_seed=args.random_seed)
    (rna_latent, atac_latent), _ = Pa.run_Pamona([X, Y])
    elapsed_time = time.time() - start_time

    print("[4/4] Saving results...")
    args.output_rna.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rna_latent, index=rna.obs_names).to_csv(args.output_rna, header=False)
    args.output_atac.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(atac_latent, index=atac.obs_names).to_csv(args.output_atac, header=False)
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": atac.shape[0] + rna.shape[0]
        }, f)


if __name__ == "__main__":
    main(parse_args())
