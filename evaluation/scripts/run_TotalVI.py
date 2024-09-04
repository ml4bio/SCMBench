r"""
Run TotalVI
"""

from scvi_colab import install
install()
import argparse
import sys, os
import pathlib
import torch
import matplotlib.pyplot as plt
import time
import numpy as np
import anndata 
import pandas as pd
import yaml
import matplotlib.pyplot as plt
import mudata as md
import muon
import scanpy as sc
import scvi


import SCMBench
from SCMBench.utils import AUTO

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
        "--input-protein", dest="input_protein", type=pathlib.Path, required=True,
        help="Path to input protein dataset (.h5ad)"
    )
    parser.add_argument(
        "-s", "--random-seed", dest="random_seed", type=int, default=0,
        help="Random seed"
    )
    parser.add_argument(
        "--train-dir", dest="train_dir", type=pathlib.Path, required=True,
        help="Path to directory where training logs and checkpoints are stored"
    )
    parser.add_argument(
        "--output-feature", dest="output_feature", type=pathlib.Path, required=True,
        help="Path of output feature latent file (.csv)"
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
    protein = anndata.read_h5ad(args.input_protein)

    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    rna.obs_names_make_unique()
    mdata = md.MuData({"rna": rna, "protein": protein})

    print("[2/4] Preprocessing...")
    start_time = time.time()
    sc.pp.highly_variable_genes(
        mdata.mod["rna"],
        n_top_genes=4000,
        flavor="seurat_v3",
        batch_key=None,
        layer="counts",
    )
    # Place subsetted counts in a new modality
    mdata.mod["rna_subset"] = mdata.mod["rna"][
        :, mdata.mod["rna"].var["highly_variable"]
    ].copy()
    mdata.update()
    scvi.model.TOTALVI.setup_mudata(
        mdata,
        rna_layer=None,
        protein_layer=None,
        batch_key=None,
        modalities={
            "rna_layer": "rna_subset",
            "protein_layer": "protein",
        },
    )

    print("[3/4] Training TotalVI...")
    vae = scvi.model.TOTALVI(mdata)
    vae.train()
    feature_latent = np.array(vae.get_latent_representation())
    elapsed_time = time.time() - start_time

    print("[4/4] Saving results...")
    args.output_feature.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(feature_latent,).to_csv(args.output_feature, header=False)
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": rna.X.shape[0]
        }, f)

if __name__ == "__main__":
    main(parse_args())