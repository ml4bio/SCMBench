#!/usr/bin/env python

r"""
Run scGPT with scVI adapter layer
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

# sys.path.append('../')
# import SCMBench
# from SCMBench.utils import AUTO

import scvi
from rich import print
import torch.nn.functional as F

# SCMBench.log.console_log_level = logging.DEBUG


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
        "--rna-pre", dest="rna_pre", type=pathlib.Path, required=True,
        help="Path to input RNA Pretrained Embeddings (.csv)"
    )
    parser.add_argument(
        "--atac-pre", dest="atac_pre", type=pathlib.Path, required=True,
        help="Path to input ATAC Pretrained Embeddings (.csv)"
    )
    parser.add_argument(
        "--output-path", dest="output_path", type=pathlib.Path, required=True,
        help="Path of output files"
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
    
    return parser.parse_args()

class PreSCVI(scvi.model.SCVI):
    def __init__(self, adata, input_key, **model_kwargs):
        super().__init__(adata, **model_kwargs)
        self.input_key = input_key

    def _get_inference_input(self, tensors):
        x = self.adata.obsm[self.input_key]
        return dict(x=x)

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors["X"]
        return dict(x=x, z=inference_outputs["qz_m"])

    def _compute_reconstruction_loss(self, generative_outputs, target, **kwargs):
        x_hat = generative_outputs["px_rate"]
        return F.mse_loss(x_hat, target)

def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    print("[1/5] Reading data...")
    rna = anndata.read_h5ad(args.input_rna)
    atac = anndata.read_h5ad(args.input_atac)
    adata = sc.concat([rna,atac],axis=1)
    adata.obs['cell_type'] = rna.obs['cell_type']

    if args.random_sleep:
        time.sleep(random.randint(0, 10))

    print("[2/5] Preprocessing...")
    start_time = time.time()
    adata.raw = adata  # keep full dimension safe

    print("[3/5] Load Pretrained Embeddings...")
    atac_emb = pd.read_csv(args.atac_pre,index_col=0, header=None)
    atac_emb=atac_emb.values
    rna_emb = pd.read_csv(args.rna_pre,index_col=0, header=None)
    rna_emb = rna_emb.values
    both=np.concatenate([rna_emb,atac_emb],axis=1)

    adata.obsm['X_FM'] = both

    print("[4/5] Training scVI...")
    PreSCVI.setup_anndata(adata)
    vae = PreSCVI(adata, input_key="X_FM", n_layers=args.n_layers, n_latent=args.n_latent, gene_likelihood="nb")
    vae.train()

    adata.obsm["X_scVI"] = vae.get_latent_representation()
    elapsed_time = time.time() - start_time

    print("[5/5] Saving results...")
    args.output_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs_names).to_csv(args.output_path / 'rna.csv', header=False)
    pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs_names).to_csv(args.output_path / 'atac.csv', header=False)
    with open(args.output_path / 'run_info.yaml','w') as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": rna.shape[0]
        }, f)

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=0.3)
    fig = sc.pl.umap(
        adata,
        color=["cell_type"],
        frameon=False,
        return_fig=True,
        show=False,
        palette='Set2',
    )
    fig.savefig(
                args.output_path /'embeddings_celltype_umap.png', dpi=300, bbox_inches="tight"
            )

if __name__ == "__main__":
    main(parse_args())
