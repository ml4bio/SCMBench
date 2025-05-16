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

import torch
import scvi
from rich import print
import scvi
from scvi.module import VAE
from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import _get_adata_minify_type
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.data.fields import (
    LayerField,
    ObsmField,
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from typing import Literal

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
        "--gene_likelihood", dest="gene_likelihood", type=str, default="zinb",
        help="gene_likelihood of the scVI model"
    )
    
    return parser.parse_args()

class CustomVAE(VAE):
    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_latent: int = 30,
        n_hidden: int = 128,
        n_layers: int = 2,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        fm_dim: int = 512,  
        dropout_rate: float = 0.1,
        dispersion: str = "gene",
        log_variational: bool = False,
        use_batch_norm: bool = True,
        bias: bool = False,
        latent_distribution: str = "normal",
        use_observed_lib_size: bool = True,
        **kwargs
    ):
        from scvi.nn import Encoder
        
        super().__init__(n_input=n_input,
                        n_batch=n_batch,
                        n_labels=n_labels,
                        n_hidden=n_hidden,
                        n_latent=n_latent,
                        dropout_rate=dropout_rate,
                        dispersion=dispersion,
                        log_variational=log_variational,
                        gene_likelihood=gene_likelihood,
                        latent_distribution=latent_distribution,
                        use_observed_lib_size=use_observed_lib_size,
                        **kwargs
            )  
        n_input_encoder = fm_dim
        self.z_encoder = Encoder(
            n_input_encoder,
            n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            distribution=latent_distribution,
            use_batch_norm=True,
            use_layer_norm=False,
            return_dist=True,
        )
        # l encoder goes from n_input-dimensional data to 1-d library size
        self.l_encoder = Encoder(
            n_input_encoder,
            1,
            n_layers=1,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=True,
            use_layer_norm=False,
            return_dist=True,
        )

    def _get_inference_input(
        self,
        tensors,
    ):
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        cont_key = REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        y_key = REGISTRY_KEYS.X_KEY
        y = tensors[y_key]

        if self.minified_data_type is None:
            x = tensors['X_FM']
            input_dict = {
                "x": x,
                "y": y,
                "batch_index": batch_index,
                "cont_covs": cont_covs,
                "cat_covs": cat_covs,
            }
        else:
            if self.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
                qzm = tensors[REGISTRY_KEYS.LATENT_QZM_KEY]
                qzv = tensors[REGISTRY_KEYS.LATENT_QZV_KEY]
                observed_lib_size = tensors[REGISTRY_KEYS.OBSERVED_LIB_SIZE]
                input_dict = {
                    "qzm": qzm,
                    "qzv": qzv,
                    "observed_lib_size": observed_lib_size,
                }
            else:
                raise NotImplementedError(f"Unknown minified-data type: {self.minified_data_type}")

        return input_dict
    
    def _regular_inference(
        self,
        x,
        y,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        n_samples=1,
    ):
        """High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = x
        if self.use_observed_lib_size:
            library = torch.log(y.sum(1)).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)
        if cont_covs is not None and self.encode_covariates:
            encoder_input = torch.cat((x_, cont_covs), dim=-1)
        else:
            encoder_input = x_
        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()
        qz, z = self.z_encoder(encoder_input, batch_index, *categorical_input)
        ql = None
        if not self.use_observed_lib_size:
            ql, library_encoded = self.l_encoder(encoder_input, batch_index, *categorical_input)
            library = library_encoded

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            if self.use_observed_lib_size:
                library = library.unsqueeze(0).expand(
                    (n_samples, library.size(0), library.size(1))
                )
            else:
                library = ql.sample((n_samples,))
        outputs = {"z": z, "qz": qz, "ql": ql, "library": library}
        return outputs
    
class PreSCVI(scvi.model.SCVI):
    _module_cls = CustomVAE
    
    def __init__(
        self,
        adata,
        **model_kwargs
    ):
        super().__init__(adata)
        
        self.module = self._module_cls(
            n_input=self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY).shape[1],
            fm_dim=self.adata_manager.get_from_registry('X_FM').shape[1],
            **model_kwargs
        )
        
        self._model_summary_string = "CustomSCVI Model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def setup_anndata(
        cls,
        adata,
        fm_key: str = "X_FM",
        layer= None,
        batch_key= None,
        labels_key= None,
        size_factor_key= None,
        categorical_covariate_keys = None,
        continuous_covariate_keys = None,
        **kwargs,
    ):
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            ObsmField(fm_key, fm_key),
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

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
    PreSCVI.setup_anndata(adata, fm_key = "X_FM")
    vae = PreSCVI(adata, n_layers=args.n_layers, n_latent=args.n_latent, gene_likelihood=args.gene_likelihood)
    vae.train()

    adata.obsm["X_scVI"] = vae.get_latent_representation()
    elapsed_time = time.time() - start_time

    print("[5/5] Saving results...")
    args.output_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs_names).to_csv(args.output_path / 'rna_embeddings.csv', header=False)
    pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs_names).to_csv(args.output_path / 'atac_embeddings.csv', header=False)
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
