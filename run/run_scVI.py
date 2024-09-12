#!/usr/bin/env python

r"""
Run scVI
"""

import argparse
import logging
import pathlib
import random
import sys
import time

import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

import SCMBench
from SCMBench.utils import AUTO

import scvi
from rich import print
from scvi.model.utils import mde

SCMBench.log.console_log_level = logging.DEBUG


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
        "--lr", dest="lr", type=float, default=2e-3,
        help="Learning rate"
    )
    parser.add_argument(
        "--random-sleep", dest="random_sleep", default=False, action="store_true",
        help="Whether to sleep random number of seconds before running, "
             "which helps distribute jobs evenly over multiple GPUs."
    )
    parser.add_argument(
        "--max-epochs", dest="max_epochs", type=int, default=300,
        help="Max training epochs"
    )
    parser.add_argument(
        "--layer", dest="layer", type=str, default="counts",
        help="Register the AnnData object with the correct key to identify the sample and the layer key with the count data"
    )
    parser.add_argument(
        "--batch_key", dest="batch_key", type=str, default="batch",
        help="Register the AnnData object with the correct key to identify the sample and the layer key with the count data"
    )
    parser.add_argument(
        "--n_layers", dest="n_layers", type=int, default=2,
        help="n_layers of the scVI model"
    )
    parser.add_argument(
        "--n_latent", dest="n_latent", type=int, default=30,
        help="n_latent of the scVI model"
    )
    parser.add_argument(
        "--gene_likelihood", dest="gene_likelihood", type=str, default="nb",
        help="gene_likelihood of the scVI model"
    )
    parser.add_argument(
        "--output-rna", dest="output_rna", type=pathlib.Path, required=True,
        help="Path of output RNA latent file (.csv)"
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
    print("[1/4] Reading data...")
    rna = anndata.read_h5ad(args.input_rna)
    if args.random_sleep:
        time.sleep(random.randint(0, 10))

    print("[2/4] Preprocessing...")
    start_time = time.time()
    rna.raw = rna  # keep full dimension safe

    print("[3/4] Training scVI...")
    scvi.model.SCVI.setup_anndata(rna)
    vae = scvi.model.SCVI(rna, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()

    rna.obsm["X_scVI"] = vae.get_latent_representation()
    elapsed_time = time.time() - start_time

    print("[4/4] Saving results...")
    args.output_rna.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rna.obsm["X_scVI"], index=rna.obs_names).to_csv(args.output_rna, header=False)
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": rna.shape[0]
        }, f)


if __name__ == "__main__":
    main(parse_args())
