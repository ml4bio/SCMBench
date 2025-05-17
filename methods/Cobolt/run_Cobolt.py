#!/usr/bin/env python

r"""
Run Cobolt
"""

import argparse
import pathlib
import sys
import time

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import yaml

import cobolt
from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt


def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-joint-rna", dest="input_joint_rna", type=pathlib.Path, required=True,
        help="Path to input joint RNA dataset (.h5ad)"
    )
    parser.add_argument(
        "--input-joint-atac", dest="input_joint_atac", type=pathlib.Path, required=True,
        help="Path to input joint ATAC dataset (.h5ad)"
    )
    parser.add_argument(  # no use
        "-s", "--random-seed", dest="random_seed", type=int, default=0,
        help="Random seed"
    )
    parser.add_argument(
        "--output-joint", dest="output_joint", type=pathlib.Path, required=True,
        help="Path of output joint modality latent file (.csv)"
    )
    parser.add_argument(
        "-r", "--run-info", dest="run_info", type=pathlib.Path, required=True,
        help="Path of output run info file (.yaml)"
    )
    parser.add_argument(
        "--unpaired", dest="unpaired", type=bool, default=False,
        help="Integrate unpaired modalities"
    )
    parser.add_argument(
        "--input-single-rna", dest="input_single_rna", type=pathlib.Path, default="",
        help="Path to input single RNA dataset (.h5ad)"
    )
    parser.add_argument(
        "--input-single-atac", dest="input_single_atac", type=pathlib.Path, default="",
        help="Path to input single ATAC dataset (.h5ad)"
    )
    parser.add_argument(
        "--output-rna", dest="output_rna", type=pathlib.Path, default="",
        help="Path of output RNA latent file (.csv)"
    )
    parser.add_argument(
        "--output-atac", dest="output_atac", type=pathlib.Path, default="",
        help="Path of output ATAC latent file (.csv)"
    )
    parser.add_argument(
        "--upper-quantile", dest="upper_quantile", type=float, default=0.99,
        help="Feature filtering parameter: upper quantile"
    )
    parser.add_argument(
        "--lower-quantile", dest="lower_quantile", type=float, default=0.7,
        help="Feature filtering parameter: lower quantile"
    )
    parser.add_argument(
        "--lr", dest="lr", type=float, default=0.001,
        help="Training parameter: learning rate"
    )
    parser.add_argument(
        "--latent-dim", dest="ld", type=int, default=10,
        help="Training parameter: latent dimension"
    )
    parser.add_argument(
        "--epoch", dest="epoch", type=int, default=20,
        help="Training parameter: epochs"
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    np.random.seed(args.random_seed)

    print("[1/4] Reading data...")

    joint_rna = anndata.read_h5ad(args.input_joint_rna)
    joint_atac = anndata.read_h5ad(args.input_joint_atac)

    if args.unpaired:
        single_rna = anndata.read_h5ad(args.input_single_rna)
        single_atac = anndata.read_h5ad(args.input_single_atac)

    # process to Cobolt required format
    joint_rna_count = joint_rna.X.tocsr().astype(float)
    joint_rna_feature = joint_rna.var_names.values.astype('str')
    joint_rna_barcode = joint_rna.obs_names.values.astype('str')
    joint_rna_cobolt = SingleData("GeneExpr", str(args.input_joint_rna).split('/')[-2], joint_rna_feature, joint_rna_count, joint_rna_barcode)

    joint_atac_count = joint_atac.X.tocsr().astype(float)
    joint_atac_feature = joint_atac.var_names.values.astype('str')
    joint_atac_barcode = joint_atac.obs_names.values.astype('str')
    joint_atac_cobolt = SingleData("ChromAccess", str(args.input_joint_atac).split('/')[-2], joint_atac_feature, joint_atac_count, joint_atac_barcode)

    if args.unpaired:
        single_rna_count = single_rna.X.tocsr().astype(float)
        single_rna_feature = single_rna.var_names.values.astype('str')
        single_rna_barcode = single_rna.obs_names.values.astype('str')
        single_rna_cobolt = SingleData("GeneExpr", str(args.input_single_rna).split('/')[-2], single_rna_feature, single_rna_count, single_rna_barcode)

        single_atac_count = single_atac.X.tocsr().astype(float)
        single_atac_feature = single_atac.var_names.values.astype('str')
        single_atac_barcode = single_atac.obs_names.values.astype('str')
        single_atac_cobolt = SingleData("ChromAccess", str(args.input_single_atac).split('/')[-2], single_atac_feature, single_atac_count, single_atac_barcode)
    
    print("[2/4] Preprocessing...")

    start_time = time.time()
    
    joint_rna_cobolt.filter_features(upper_quantile=args.upper_quantile, lower_quantile=args.lower_quantile)
    joint_atac_cobolt.filter_features(upper_quantile=args.upper_quantile, lower_quantile=args.lower_quantile)
    if (args.unpaired):
        single_rna_cobolt.filter_features(upper_quantile=args.upper_quantile, lower_quantile=args.lower_quantile)
        single_atac_cobolt.filter_features(upper_quantile=args.upper_quantile, lower_quantile=args.lower_quantile)

    if args.unpaired:
        multi_dt = MultiomicDataset.from_singledata(joint_rna_cobolt, joint_atac_cobolt, single_rna_cobolt, single_atac_cobolt)
    else:
        multi_dt = MultiomicDataset.from_singledata(joint_rna_cobolt, joint_atac_cobolt)
    print(multi_dt)
    
    print("[3/4] Training Cobolt...")

    model = Cobolt(dataset=multi_dt, lr=args.lr, n_latent=args.ld)
    model.train(num_epochs=args.epoch)
    model.calc_all_latent()

    latent, _ = model.get_all_latent()
    if args.unpaired:
        rna_latent = latent[2*joint_rna.shape[0]:2*joint_rna.shape[0]+single_rna.shape[0]]  # embeddings for RNAseq cells
        atac_latent = latent[3*joint_rna.shape[0]+single_rna.shape[0]:]  # embeddings for ATACseq cells
    else:
        joint_latent = latent

    elapsed_time = time.time() - start_time
    print(joint_latent.shape)
    
    print("[4/4] Saving results...")

    if args.unpaired:
        args.output_rna.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rna_latent, index=single_rna.obs_names).to_csv(args.output_rna, header=False)
        args.output_atac.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(atac_latent, index=single_atac.obs_names).to_csv(args.output_atac, header=False)
    else:
        args.output_joint.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(joint_latent, index=joint_rna.obs_names).to_csv(args.output_joint, header=False)
    
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    n_cells = joint_rna.shape[0]
    if args.unpaired:
        n_cells += single_rna.shape[0] + single_atac.shape[0]
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": n_cells
        }, f)


if __name__ == "__main__":
    main(parse_args())
