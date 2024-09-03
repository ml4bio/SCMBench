#!/usr/bin/env python

r"""
Run scJoint
"""

import argparse
import pathlib
import sys
import time
import os

import anndata
import pandas as pd
import scipy.sparse
import yaml
import numpy as np
import torch

sys.path.append("../../../custom")
from scJoint.config import Config
from scJoint.process_db import label_parsing
from scJoint.trainingprocess_stage1 import TrainingProcessStage1
from scJoint.knn import KNN
from scJoint.trainingprocess_stage3 import TrainingProcessStage3


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
        help="Path to input ATAC dataset in gene activity matrix format (.csv)"
    )
    parser.add_argument(  # no use
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
    parser.add_argument(
        "--lr-stage1", dest="lr_stage1", type=float, default=0.01,
        help="Training parameter: learning rate for stage 1"
    )
    parser.add_argument(
        "--lr-stage3", dest="lr_stage3", type=float, default=0.01,
        help="Training parameter: learning rate for stage 3"
    )
    parser.add_argument(
        "--lr-decay-epoch", dest="lr_decay_epoch", type=int, default=20,
        help="Training parameter: number of epochs learning rate decay"
    )
    parser.add_argument(
        "--epoch-stage1", dest="epoch_stage1", type=int, default=20,
        help="Training parameter: number of epochs for stage 1"
    )
    parser.add_argument(
        "--epoch-stage3", dest="epoch_stage3", type=int, default=20,
        help="Training parameter: number of epochs for stage 3"
    )
    parser.add_argument(
        "--latent-dim", dest="ld", type=int, default=64,
        help="Training parameter: latent dimension"
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    torch.set_num_threads(1)
    os.environ['OMP_NUM_THREADS'] = '1'
    np.random.seed(args.random_seed)
    torch.manual_seed(args.random_seed)
    
    print("[1/4] Reading data...")

    rna = anndata.read_h5ad(args.input_rna)
    atac = pd.read_csv(args.input_atac, index_col=0, float_precision='round_trip') # gene activity matrix
    atac_obs = anndata.read_h5ad(str(args.input_atac).replace('ATAC_GAM.csv', 'RNA.h5ad')).obs
    atac = anndata.AnnData(X=scipy.sparse.csr_matrix(atac.values), obs=atac_obs, var=pd.DataFrame(index=atac.columns), 
                           dtype=scipy.sparse.csr_matrix(atac.values).dtype)

    # process to scJoint required format
    common_genes = rna.var_names.intersection(atac.var_names)
    rna = rna[:, common_genes]
    atac = atac[:, common_genes]
    print('Number of common genes:', rna.shape[1])
    if rna.shape[1] == 0:
        print('TOO FEW GENES!')
        exit()

    common_cell_types = set(rna.obs['cell_type']).intersection(set(atac.obs['cell_type']))
    rna = rna[rna.obs['cell_type'].isin(common_cell_types)]
    atac = atac[atac.obs['cell_type'].isin(common_cell_types)]
    print('Number of cells of common types:', rna.shape[0])
    if rna.shape[0] == 0:
        print('TOO FEW CELLS!')
        exit()

    config = Config(args, rna.shape[1], len(common_cell_types))

    rna_sparse = rna.X.tocsr()
    scipy.sparse.save_npz(config.rna_paths[0], rna_sparse)
    atac_sparse = atac.X.tocsr() #
    scipy.sparse.save_npz(config.atac_paths[0], atac_sparse) #
    rna.obs['cell_type'].to_csv(config.rna_labels[0].replace('txt','csv'))
    label_parsing([config.rna_labels[0].replace('txt','csv')], []) #

    print("[2/4] Preprocessing...")

    start_time = time.time()
    
    print("[3/4] Training scJoint...")

    # stage1 training
    print('Training start [Stage1]')
    model_stage1= TrainingProcessStage1(config)    
    for epoch in range(config.epochs_stage1):
        print('Epoch:', epoch)
        model_stage1.train(epoch)
    print('Writing embeddings and predictions in Stage 1')
    model_stage1.write_embeddings()
    print('Stage 1 finished')
    
    # KNN
    print('KNN [Stage2]')
    KNN(config, neighbors = 30, knn_rna_samples=20000)
    print('Stage 2 finished')
    
    # stage3 training
    print('Training start [Stage3]')
    model_stage3 = TrainingProcessStage3(config)    
    for epoch in range(config.epochs_stage3):
       print('Epoch:', epoch)
       model_stage3.train(epoch)
    print('Stage 3 finished')

    elapsed_time = time.time() - start_time
    
    print("[4/4] Saving results...")

    print('Write embeddings [Stage3]')
    model_stage3.write_embeddings(args.output_rna, args.output_atac)
    
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": rna.shape[0] + atac.shape[0]
        }, f)


if __name__ == "__main__":
    main(parse_args())
