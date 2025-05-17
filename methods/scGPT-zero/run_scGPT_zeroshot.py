r"""
Run scGPT zeroshot
"""
import sys
import warnings
import argparse
import pathlib
from pathlib import Path
import os 

import scanpy as sc
import pandas as pd

sys.path.insert(0, "../")

import scgpt as scg
import matplotlib.pyplot as plt

def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-path", dest="input_path", type=pathlib.Path, required=True,
        help="Path to input RNA or ACTIVATE(ATAC) dataset (.h5ad)"
    )
    parser.add_argument(
        "--model-path", dest="model_path", type=pathlib.Path, required=True,
        help="Path to scGPT pretrained model directory path"
    )
    parser.add_argument(
        "--output-path", dest="output_path", type=pathlib.Path, required=True,
        help="Path of output (.csv)"
    )    

    return parser.parse_args()

def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    print("[1/3] Reading data...")
    plt.style.context('default')
    warnings.simplefilter("ignore", ResourceWarning)

    model_dir = args.model_path
    smaple_data_path=args.input_path

    adata = sc.read_h5ad(smaple_data_path)
    adata.var['gene_name']=list(adata.var.index)

    gene_col = "gene_name"
    cell_type_key = "cell_type"

    celltype_id_labels = adata.obs[cell_type_key].astype("category").cat.codes.values
    adata = adata[celltype_id_labels >= 0]

    print("[2/3] Encoding data...")
    embed_adata = scg.tasks.embed_data(
        adata,
        model_dir,
        gene_col=gene_col,
        batch_size=64,
    )

    print("[3/3] Saving embeddings...")
    emb=pd.DataFrame(embed_adata.obsm['X_scGPT'],index=adata.obs_names)
    os.makedirs(args.output_path.parent,exist_ok=True)
    emb.to_csv(args.output_path,header=False)

if __name__ == "__main__":
    main(parse_args())