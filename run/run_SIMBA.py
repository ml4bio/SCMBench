#!/usr/bin/env python

r"""
Run SIMBA
"""

import argparse
import os
import pathlib
import sys
import time
from typing import Optional, Union

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import sklearn.preprocessing
import sklearn.utils.extmath
import yaml
import simba as si
import torch
from signal import signal, SIGPIPE, SIG_DFL  
signal(SIGPIPE,SIG_DFL)


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
    #parser.add_argument(
    #    "-g", "--graph-save-dir", dest="graph_save_dir", type=pathlib.Path, required=True,
    #    help="Path of graph stats and pbg config"
    #)
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    np.random.seed(args.random_seed)

    print("[1/4] Reading data...")
    rna = anndata.read_h5ad(args.input_rna)
    atac = anndata.read_h5ad(args.input_atac)

    # rename to simba required format
    atac.var['chr'] = atac.var['chrom']
    atac.var['start'] = atac.var['chromStart']
    atac.var['end'] = atac.var['chromEnd']

    print("[2/4] Preprocessing...")
    start_time = time.time()
    
    si.pp.filter_peaks(atac, min_n_cells=3)
    si.pp.cal_qc_atac(atac)
    si.pp.pca(atac, n_components=50)
    si.pp.select_pcs_features(atac)

    si.pp.filter_genes(rna, min_n_cells=3)
    si.pp.cal_qc_rna(rna)
    si.pp.normalize(rna,method='lib_size')
    si.pp.log_transform(rna)
    si.pp.select_variable_genes(rna, n_top_genes=4000)
    # discretize RNA expression
    si.tl.discretize(rna, n_bins=5)

    # Infer edges between cells of different modalities
    adata_rna_atac = si.tl.gene_scores(atac, genome='hg38', use_gene_weigt=True, use_top_pcs=True)
    si.pp.filter_genes(adata_rna_atac, min_n_cells=3)
    si.pp.cal_qc_rna(adata_rna_atac) 
    si.pp.normalize(adata_rna_atac, method='lib_size')
    si.pp.log_transform(adata_rna_atac)

    adata_CrnaCatac = si.tl.infer_edges(rna, adata_rna_atac, n_components=15, k=15)

    
    print("[3/4] Training SIMBA...")

    si.tl.gen_graph(list_CP=[atac],
                list_CG=[rna],
                list_CC=[adata_CrnaCatac],
                copy=False,
                use_highly_variable=True,
                use_top_pcs=True,
                dirname='graph0')
    dict_config = si.settings.pbg_params.copy()
    if torch.cuda.is_available():
        dict_config['num_gpus'] = 1
        dict_config['workers'] = 0
    else:
        dict_config['num_gpus'] = 0
    print("settings:", dict_config)
    si.tl.pbg_train(pbg_params=dict_config, auto_wd=True, save_wd=True, output='model')

    si.load_graph_stats()
    si.load_pbg_config()
    dict_adata = si.read_embedding()

    atac_latent = dict_adata['C']  # embeddings for ATACseq cells
    rna_latent = dict_adata['C2']  # embeddings for RNAseq cells
    #adata_G = dict_adata['G']  # embeddings for genes
    #adata_P = dict_adata['P']  # embeddings for peaks
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
