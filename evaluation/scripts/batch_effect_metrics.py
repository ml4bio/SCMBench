#!/usr/env/bin python

r"""
Compute cell integration metrics
"""

import argparse
import functools
import pathlib
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import yaml

from scipy.sparse import issparse

import SCMBench
import SCMBench.metrics
import SCMBench.batch_correction_metrics

def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Compute integration metrics for paired samples"
    )
    parser.add_argument(
        "-d", "--datasets", dest="datasets", type=pathlib.Path, required=True,
        nargs="+", help="Path to datasets (.h5ad)"
    )
    parser.add_argument(
        "-l", "--latents", dest="latents", type=pathlib.Path, required=True,
        nargs="+", help="Path to latent embeddings (.csv)"
    )
    parser.add_argument(
        "--cell-type", dest="cell_type", type=str, default="cell_type",
        help="Column name in obs specifying cell types"
    )
    parser.add_argument(
        "--domain", dest="domain", type=str, default="domain",
        help="Column name in obs specifying domain"
    )
    parser.add_argument(
        "-p", "--paired", dest="paired", default=False, action="store_true",
        help="Whether the latent embeddings are paired"
    )
    parser.add_argument(
        "-o", "--output", dest="output", type=pathlib.Path, required=True,
        help="Path to output file (.yaml)"
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    if len(args.datasets) != len(args.latents):
        raise RuntimeError("Datasets and latents should have the same number of entries!")

    print("[1/3] Reading data...")
    datasets = [anndata.read_h5ad(item) for item in args.datasets]
    latents = [pd.read_csv(item, header=None, index_col=0).to_numpy() for item in args.latents]
    for i in range(len(datasets)):
        datasets[i].obsm['X_embed'] = latents[i]
    print("[2/3] Computing metrics...")
    masks = [np.apply_along_axis(lambda x: ~np.any(np.isnan(x)), 1, latent) for latent in latents]
    for i, mask in enumerate(masks):
        rm_pct = 100 * (1 - mask.sum() / mask.size)
        if rm_pct:
            print(f"Ignoring {rm_pct:.1f}% cells in dataset {i} due to missing values!")
    
    combined_adata = anndata.concat([dataset[mask] for dataset, mask in zip(datasets, masks)])
    combined_adata.obs['batch'] = combined_adata.obs['batch'].astype(str).astype('category')
    combined_adata.obs['cell_type'] = combined_adata.obs['cell_type'].astype(str).astype('category')
    sc.pp.neighbors(combined_adata, use_rep="X_embed")

    metrics = {
        "avg_silhouette_width_batch":
            SCMBench.batch_correction_metrics.silhouette_batch(combined_adata, batch_key="batch", label_key="cell_type", embed="X_embed"),
        "ilisi":
            SCMBench.batch_correction_metrics.ilisi_graph(combined_adata, batch_key="batch", type_="embed", use_rep="X_embed"),
        "kbet":
            SCMBench.batch_correction_metrics.kBET(combined_adata, batch_key="batch", label_key="cell_type", type_="embed", embed="X_embed"),
        "graph_connectivity":
            SCMBench.batch_correction_metrics.graph_connectivity(combined_adata, label_key="cell_type"),
    }
    
    # round results to .4f
    for k, v in metrics.items():
        metrics[k] = round(float(v), 4)
    print("[3/3] Saving results...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        yaml.dump(metrics, f)


if __name__ == "__main__":
    main(parse_args())