from mofapy2.run.entry_point import entry_point
import h5py

import argparse
import time

import argparse
import pathlib
import sys
import time
import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
sys.path.append("../../..")
sys.path.append("../custom")
import SCMBench

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
    ent = entry_point()

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

    print("[3/4] Training MOFA+...")
    output_dir = args.run_info.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    X = rna.obsm["X_pca"]
    Y = atac.obsm["X_lsi"]

    data_tr1 = X
    data_tr2 = Y


    D = [data_tr1.shape[1], data_tr2.shape[1]] # Number of features per view
    M = len(D)      # Number of views
    K = 5           # Number of factors
    N = [data_tr1.shape[0]]   # Number of samples per group
    G = len(N)


    data_mat = [[None for g in range(G)] for m in range(M)]

    data_mat[0][0] = data_tr1
    data_mat[1][0] = data_tr2
            

    ent.set_data_options(
        scale_groups = False, 
        scale_views = False
    )


    ent.set_data_matrix(data_mat, likelihoods = ["gaussian","gaussian"])


    ent.set_model_options(
        factors = 100, 
        spikeslab_weights = True, 
        ard_factors = True,
        ard_weights = True
    )


    ent.set_train_options(
        iter = 100, 
        convergence_mode = "fast", 
        startELBO = 1, 
        freqELBO = 1, 
        dropR2 = 0.001, 
        gpu_mode = True, 
        verbose = False, 
        seed = 1
    )


    ent.build()
    ent.run()
    if os.path.exists('./MOFAmodels/tmp.hdf5'):
        os.remove('./MOFAmodels/tmp.hdf5')

    ent.save(outfile='./MOFAmodels/tmp.hdf5')

    #end-start
    f = h5py.File('./MOFAmodels/tmp.hdf5', "r")

    rna_latent=np.array(f['data']['view0']['group0'])

    atac_latent=np.array(f['data']['view1']['group0'])

    rna_W=np.array(f['expectations']['W']['view0'])

    atac_W=np.array(f['expectations']['W']['view1'])

    rna_latent=np.dot(rna_latent,rna_W.transpose())
    atac_latent=np.dot(atac_latent,atac_W.transpose())
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