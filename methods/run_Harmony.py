#!/usr/bin/env python

import harmonypy as hm
import scanpy as sc
import pandas as pd
import numpy as np
import time
import yaml
import argparse
import pathlib
import logging
import pathlib
import random
import sys
import time


# args = argparse.Namespace(
#     input_rna = '/ailab/user/chenpengan/zgy/DualOmics/Dual_datasets/Muto-2021-small/Muto-2021-small-RNA.h5ad',
#     input_atac = '/ailab/user/chenpengan/zgy/DualOmics/Dual_datasets/Muto-2021-small/Muto-2021-small-ACTIVE.h5ad',
#     random_seed = 666,
#     output_rna = "/ailab/user/chenpengan/zgy/DualOmics/harmony/Muto-2021-small-RNA_latent.csv",
#     output_atac = "/ailab/user/chenpengan/zgy/DualOmics/harmony/Muto-2021-small-ACTIVE_latent.csv"
# )

def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s","--random_seed", dest="random_seed", type=int, default = 0,
        help="Random seed"
    )
    parser.add_argument(
        "--input-rna", dest="input_rna", type=pathlib.Path, required=True,
        help="Path to input RNA dataset (.h5ad)"
    )
    parser.add_argument(
        "--input-atac", dest="input_atac", type=pathlib.Path, required=True,
        help="Path to input ATAC dataset converted to RNA (.h5ad)"
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
    # Set random seed for reproducibility
    np.random.seed(args.random_seed)

    print("[1/4] Reading data...")

    # Load the RNA and ATAC data
    rna = sc.read(args.input_rna)
    atac = sc.read(args.input_atac)
    print('raw:',rna.shape, atac.shape)


        # Ensure the genes (var) match between RNA and ATAC
        # assert np.all(atac.var_names == rna.var_names),"Gene names do not match between RNA and ATAC data!"

    # Adjust row names to avoid collisions
    rna.obs_names = [f"{name}.RNA" for name in rna.obs_names]
    atac.obs_names = [f"{name}.ATAC" for name in atac.obs_names]

    # Extract highly variable genes
    hvg = rna.var_names[rna.var['highly_variable']]

    # Combine both datasets
    combined = rna.concatenate(atac, join='outer', batch_key='domain')

    # Preprocessing: Normalize data and scale
    print("[2/4] Data preprocessing...")
    start_time = time.time()
    sc.pp.normalize_total(combined)
    sc.pp.log1p(combined)
    sc.pp.scale(combined)
    combined.raw = combined  # Preserve raw counts
    sc.tl.pca(combined, random_state=args.random_seed, n_comps=50)

    # Run Harmony
    print("[3/4] Running Harmony...")
    vars_use = ['domain']

    ho  = hm.run_harmony(combined.obsm['X_pca'], combined.obs,vars_use)

    print('combined',combined.shape)
    res = pd.DataFrame(ho.Z_corr).T
    res.index = combined.obs.index

    elapsed_time = time.time() - start_time

    print("[4/4] Saving results...")

    # Save results: Extract the embeddings from Harmony

    rna_latent = res[res.index.str.endswith('.RNA-0')]
    atac_latent = res[res.index.str.endswith('.ATAC-1')]
    

    # Remove suffix from row names
    rna_latent.index = [name.replace(".RNA-0", "") for name in rna_latent.index]
    atac_latent.index = [name.replace(".ATAC-1", "") for name in atac_latent.index]

    print(rna_latent.shape,atac_latent.shape)
    # Save to CSV
    rna_latent.to_csv(args.output_rna, header=False)
    atac_latent.to_csv(args.output_atac, header=False)

    # Save the runtime and cell count information
    run_info = {
        "args": vars(args),  # Convert argparse Namespace to dict
        "time": elapsed_time,
        "n_cells": rna.n_obs + atac.n_obs
    }

    with open(args.run_info, 'w') as f:
        yaml.dump(run_info, f)

if __name__ == "__main__":
    main(parse_args())
