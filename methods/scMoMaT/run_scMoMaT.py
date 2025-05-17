#!/usr/bin/env python

r"""
Run scMoMaT
"""
import argparse
import sys, os
import pathlib
import torch
import matplotlib.pyplot as plt
import scipy.sparse as sp
import logging
import random
import time
import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

import SCMBench
from SCMBench.utils import AUTO

import scmomat 

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
        "--latent-dimensions", dest="latent_dimensions", type=int, default=30,
        help="Number of latent dimensions, key hyper-parameter, 20~30 works for most of the cases."
    )
    parser.add_argument(
        "--lamb", dest="lamb", type=float, default=0.001,
        help="Weight on regularization term, default value."
    )
    parser.add_argument(
        "--T", dest="num_iteration", type=int, default=4000,
        help="Number of total iterations, default value."
    )
    parser.add_argument(
        "--TT", dest="retraining_iteration", type=int, default=2000,
        help="Number of total iterations for retraining, default value."
    )
    parser.add_argument(
        "--interval", dest="interval", type=int, default=1000,
        help="Print the result after each ``interval'' iterations, default value."
    )
    parser.add_argument(
        "--batch-size", dest="batch_size", type=float, default=0.1,
        help="Batch size for each iteraction, default value."
    )
    parser.add_argument(
        "--lr", dest="lr", type=float, default=1e-2,
        help="Learning rate, default value."
    )
    parser.add_argument(
        "--n-batches", dest="n_batches", type=int, default=1,
        help="Number of batches."
    )
    parser.add_argument(
        "--train-dir", dest="train_dir", type=pathlib.Path, required=True,
        help="Path to directory where training logs and checkpoints are stored"
    )
    parser.add_argument(
        "--output-marker-feature", dest="output_marker_feature_score", type=pathlib.Path, required=True,
        help="Path of output marker feature score file (.csv)"
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
    atac = anndata.read_h5ad(args.input_atac)
    

    print("[2/4] Preprocessing...")
    start_time = time.time()
    counts_rna = np.array(rna.X.copy().todense())
    counts_rna = scmomat.preprocess(counts_rna, modality = "RNA", log = False)
    counts_atac = np.array(atac.X.copy().todense())
    counts_atac = scmomat.preprocess(counts_atac, modality = "ATAC")

    # obtain the feature name
    genes = rna.var_names
    regions = atac.var_names
    feats_name = {"rna": genes, "atac": regions}
    counts_rnas = [counts_rna]
    counts_atacs = [counts_atac]

    # CREATE THE COUNTS OBJECT
    counts = {"feats_name": feats_name, "nbatches": args.n_batches , "rna":counts_rnas, "atac": counts_atacs}

    print("[3/4] Training scMoMaT...")

    # 1st stage training, learning cell factors
    model = scmomat.scmomat_model(counts = counts, K = args.latent_dimensions, interval = args.interval, lr = args.lr, lamb = args.lamb, seed = args.random_seed, batch_size = args.batch_size)
    losses = model.train_func(T = args.num_iteration)
    # extract cell factors
    feature_latent = model.extract_cell_factors()
    feature_latent_save = np.squeeze(feature_latent.copy())
    n_neighbors = 100
    r = None
    knn_indices, knn_dists = scmomat.calc_post_graph(feature_latent, n_neighbors, njobs = 8, r = r)
    resolution = 0.9
    labels_leiden = scmomat.leiden_cluster(X = None, knn_indices = knn_indices, knn_dists = knn_dists, resolution = resolution)

    # Retraining
    T = args.retraining_iteration
    model2 = scmomat.scmomat_retrain(model = model, counts =  counts, labels = labels_leiden, lamb = args.lamb)
    losses = model2.train(T = T)
    elapsed_time = time.time() - start_time

    # extract marker feature scores
    C_feats = model2.extract_marker_scores()
    C_gene = C_feats["rna"]

    print("[4/4] Saving results...")
    args.output_marker_feature_score.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(C_gene,).to_csv(args.output_marker_feature_score, header=False)
    args.output_feature.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(feature_latent_save,).to_csv(args.output_feature, header=False)
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": atac.shape[0]
        }, f)

        

if __name__ == "__main__":
    main(parse_args())