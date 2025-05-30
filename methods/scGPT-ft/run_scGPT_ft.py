r"""
Run scGPT finetuning with pretrained transformer module 
"""

import copy
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

from pathlib import Path
import sys
import time
import warnings
import pathlib

import torch
from anndata import AnnData
import anndata
import scanpy as sc
import numpy as np
import wandb
from scipy.sparse import issparse
from torch import nn
from sklearn.model_selection import train_test_split
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)
import pandas as pd
import scipy.sparse as sp
from sklearn import preprocessing
import yaml 
import time

sys.path.append('../methods/source_code/')
sys.path.insert(0, "../")
from scgpt import prepare_data, prepare_dataloader, define_wandb_metrcis, evaluate, eval_testdata, train
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.model import MultiOmicTransformerModel

import scgpt as scg
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.tokenizer import random_mask_value
from scgpt.loss import (
    masked_mse_loss,
    masked_relative_error,
    criterion_neg_log_bernoulli,
)
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed, category_str2int, eval_scib_metrics
from scgpt_model import eval_testdata_save_embed
import argparse

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
        "--model-path", dest="model_path", type=pathlib.Path, required=True,
        help="Path to scGPT pretrained model directory path"
    )
    parser.add_argument(
        "--output-dir", dest="output_dir", type=pathlib.Path, required=True,
        help="Directory of output files"
    )    
    parser.add_argument(
        "--pretrain", dest="pretrain", type=bool, default=False, 
        help="if with (True) or without (False) pretrained transformer"
    )
    parser.add_argument(
        "--load_layers", dest="load_layers", type=int, default=0, 
        help="number of transformer layers loaded to finetune"
    )
    parser.add_argument(
        "--nlayers", dest="nlayers", type=int, default=4, 
        help="number of transformer layers to initialize"
    )
    return parser.parse_args()

def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """
    sc.set_figure_params(figsize=(4, 4))
    os.environ["KMP_WARNINGS"] = "off"
    warnings.filterwarnings('ignore')

    print("[1/6] Setting Config...")
    hyperparameter_defaults = dict(
        task = 'multiomic',
        seed=42,
        dataset_name='SCMBench', # Dataset name "BMMC"
        do_train=True, # Flag to indicate whether to do update model parameters during training
        load_model=args.model_path, # Path to pre-trained model
        # freeze = True, #freeze
        # freeze = False,
        load_layers = args.load_layers,
        GEP=True, # Gene expression modelling
        GEPC=True, # Gene expression modelling for cell objective
        CLS=False,
        ESC=False,
        DAR = True, # DAR objective weight for batch correction
        DSBN = False,  # Domain-spec batchnorm,
        mask_ratio=0.4, # Default mask ratio
        explicit_zero_prob = False,  # whether explicit bernoulli for zeros
        ecs_thres=0,  # Elastic cell similarity objective, 0.0 to 1.0, 0.0 to disable
        dab_weight=1.0,
        use_batch_labels = True,
        use_mod = True,
        per_seq_batch_sample = False,
        epochs=25, # Default number of epochs for fine-tuning
        input_layer_key = "X_binned", # Default expression value binning in data pre-processing
        n_bins=51, # Default number of bins for value binning in data pre-processing
        n_hvg = 1200,  # Default number of highly variable genes
        n_hvp = 4000,
        max_seq_len = 4001, # # Default n_hvg+1
        lr=1e-3, # Default learning rate for fine-tuning
        batch_size=16, # Default batch size for fine-tuning
        layer_size=512,
        nlayers=args.nlayers,
        nhead=8, # if load model, batch_size, layer_size, nlayers, nhead will be ignored
        dropout=0.2, # Default dropout rate during model fine-tuning
        schedule_ratio=0.95,  # Default rate for learning rate decay
        save_eval_interval=5, # Default model evaluation interval
        log_interval=100, # Default log interval
        fast_transformer=False, # Default setting
        pre_norm=False, # Default setting
        amp=True,  # Default setting: Automatic Mixed Precision
        pad_token = "<pad>",
        mask_value = -1,
        pad_value = -2,
        include_zero_gene = False,
    )

    run = wandb.init(
        config=hyperparameter_defaults,
        project="scGPT",
        reinit=True,
        settings=wandb.Settings(start_method="fork"),
        mode="offline"
    )
    config = wandb.config
    print(config)

    set_seed(config.seed)

    print("[2/6] Reading data...")
    special_tokens = [config.pad_token, "<cls>", "<eoc>"]
    dataset_name = config.dataset_name
    save_dir = args.output_dir
    save_dir.mkdir(parents=True, exist_ok=True)
    print(f"save to {save_dir}")
    logger = scg.logger
    scg.utils.add_file_handler(logger, save_dir / "run.log")

    input_rna_dir=args.input_rna
    input_atac_dir=args.input_atac
    output_rna_dir=str(save_dir)+'/'+'rna_embeddings.csv'
    output_atac_dir=str(save_dir)+'/'+'atac_embeddings.csv'
    
    if dataset_name == 'SCMBench':
        adata = anndata.read_h5ad(input_rna_dir)
        adata_protein = anndata.read_h5ad(input_atac_dir)
        adata.obs["celltype"] = adata.obs["cell_type"].astype(str).astype('category')
        adata.var["gene_name"] = adata.var.index.tolist()
        le = preprocessing.LabelEncoder()
        encoded_batch = le.fit_transform(adata.obs['domain'].values)
        adata.obs["batch_id"] =  encoded_batch
        adata.obs["str_batch"] = adata.obs["batch_id"].astype('category')
        data_is_raw = True
        print('#')

    if dataset_name == 'SCMBench-simi':
        adata = anndata.read_h5ad(input_rna_dir)
        adata_protein = anndata.read_h5ad(input_atac_dir)
        adata.obs["celltype"] = adata.obs["cell_type"].astype(str).astype('category')
        adata.var["gene_name"] = adata.var.index.tolist()
        le = preprocessing.LabelEncoder()
        encoded_batch = le.fit_transform(adata.obs['batch'].values)
        adata.obs["batch_id"] =  encoded_batch
        adata.obs["str_batch"] = adata.obs["batch_id"].astype('category')
        data_is_raw = True
        print('#')

    if dataset_name == 'BMMC':
        adata = sc.read('scmdata/BMMC_processed.h5ad')
        # subset to first 3 donors with B, Mono and T cell subtypes
        adata = adata[adata.obs.DonorID.isin([10886, 11466, 12710]) & adata.obs.cell_type.isin(np.unique(adata.obs.cell_type.values)[:17])]
        adata.obs["celltype"] = adata.obs["cell_type"].astype(str).astype('category')
        adata.var["gene_name"] = adata.var.index.tolist()
        le = preprocessing.LabelEncoder()
        encoded_batch = le.fit_transform(adata.obs['batch'].values)
        adata.obs["batch_id"] =  encoded_batch
        adata.obs["str_batch"] = adata.obs["batch_id"].astype('category')
        adata_protein = adata[:, adata.var.feature_types.isin(['ADT'])].copy()
        adata_protein.var.index = ['p_' + i for i in adata_protein.var.index]
        adata = adata[:, adata.var.feature_types.isin(['GEX'])].copy()
        data_is_raw = False

    print("[3/6] Preprocessing...")
    start_time = time.time()
    if config.use_mod:
        gene_rna_df = pd.DataFrame(index = adata.var.index.tolist())
        gene_rna_df['mod'] = 'RNA'
        gene_protein_df = pd.DataFrame(index = adata_protein.var.index.tolist())
        gene_protein_df['mod'] = 'Protein'
        gene_loc_df = pd.concat([gene_rna_df, gene_protein_df])
        gene_loc_df['mod'] = gene_loc_df['mod'].astype('category')

    if config.load_model is not None:
        model_dir = Path(config.load_model)
        model_config_file = model_dir / "args.json"
        model_file = model_dir / "best_model.pt"
        vocab_file = model_dir / "vocab.json"

        vocab = GeneVocab.from_file(vocab_file)
        for s in special_tokens:
            if s not in vocab:
                vocab.append_token(s)

        adata.var["id_in_vocab"] = [
            1 if gene in vocab else -1 for gene in adata.var["gene_name"]
        ]
        gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
        logger.info(
            f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
            f"in vocabulary of size {len(vocab)}."
        )
        old_vocab = vocab
    embsize = config.layer_size
    nhead = config.nhead
    nlayers = config.nlayers
    d_hid = config.layer_size

    preprocessor = Preprocessor(
        use_key="X",  # the key in adata.layers to use as raw data
        filter_gene_by_counts=1,  # step 1
        filter_cell_by_counts=1,  # step 2
        normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
        result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
        log1p=data_is_raw,  # 4. whether to log1p the normalized data
        result_log1p_key="X_log1p",
        subset_hvg=config.n_hvg,  # 5. whether to subset the raw data to highly variable genes
        hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
        binning=config.n_bins,  # 6. whether to bin the raw data and to what number of bins
        result_binned_key="X_binned",  # the key in adata.layers to store the binned data
    )
    preprocessor(adata, batch_key=None)

    preprocessor_protein = Preprocessor(
        use_key="X",  # the key in adata.layers to use as raw data
        filter_gene_by_counts=0,  # step 1
        filter_cell_by_counts=False,  # step 2
        normalize_total=False,  # 3. whether to normalize the raw data and to what sum
        result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
        log1p=data_is_raw,  # 4. whether to log1p the normalized data
        result_log1p_key="X_log1p",
        subset_hvg=False,  # 5. whether to subset the raw data to highly variable genes
        hvg_flavor=None,
        binning=config.n_bins,  # 6. whether to bin the raw data and to what number of bins
        result_binned_key="X_binned",  # the key in adata.layers to store the binned data
    )
    preprocessor_protein(adata_protein, batch_key=None)

    data_combined = np.concatenate([adata.layers["X_binned"], adata_protein.layers["X_binned"]], axis=1)
    adata = AnnData(
        X=data_combined,
        obs=adata.obs,
        var=pd.DataFrame(index=adata.var_names.tolist() + adata_protein.var_names.tolist()),
        layers={"X_binned": data_combined,}
    )
    adata.var["gene_name"] = adata.var.index.tolist()

    if config.per_seq_batch_sample:
        # sort the adata by batch_id in advance
        adata_sorted = adata[adata.obs["batch_id"].argsort()].copy()

    all_counts = (
        adata.layers[config.input_layer_key].A
        if issparse(adata.layers[config.input_layer_key])
        else adata.layers[config.input_layer_key]
    )
    genes = adata.var["gene_name"].tolist()

    celltypes_labels = adata.obs["cell_type"].tolist()  # make sure count from 0
    num_types = len(set(celltypes_labels))
    celltypes_labels = np.array(celltypes_labels)

    batch_ids = adata.obs["batch_id"].tolist()
    num_batch_types = len(set(batch_ids))
    batch_ids = np.array(batch_ids)

    if config.use_mod:
        mod_type = np.array([gene_loc_df.loc[g, 'mod'] for g in genes])
        vocab_mod = Vocab(VocabPybind(np.unique(gene_loc_df['mod']).tolist() + special_tokens, None))
        vocab_mod.set_default_index(vocab_mod["<pad>"])
        mod_type = np.array(vocab_mod(list(mod_type)), dtype=int)
        ntokens_mod = len(vocab_mod)

    (
        train_data,
        valid_data,
        train_celltype_labels,
        valid_celltype_labels,
        train_batch_labels,
        valid_batch_labels,
    ) = train_test_split(
        all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
    )


    num_of_non_zero_genes = [
        np.count_nonzero(train_data[i]) for i in range(train_data.shape[0])
    ]
    print(f"max num of non_zero genes: {np.max(num_of_non_zero_genes)}")
    print(f"min num of non_zero genes: {np.min(num_of_non_zero_genes)}")
    print(f"average num of non_zero genes: {np.mean(num_of_non_zero_genes)}")
    print(
        f"99% quantile num of non_zero genes: {np.quantile(num_of_non_zero_genes, 0.99)}"
    )
    print(f"max original values: {np.max(train_data)}")
    print(
        f"average original non_zero values: {np.mean(train_data[np.nonzero(train_data)])}"
    )
    print(
        f"99% quantile original non_zero values: {np.quantile(train_data[np.nonzero(train_data)], 0.99)}"
    )
    print(f"num of celltypes: {num_types}")

    print("[4/6] Loading model and tokenizing...")
    if config.load_model is None:
        vocab = Vocab(VocabPybind(genes + special_tokens, None))
        vocab.set_default_index(vocab["<pad>"])
        gene_ids = np.array(vocab(genes), dtype=int)
    else:
        pretrained_genes = [g for g in genes + special_tokens if g in old_vocab]
        new_genes = [g for g in genes + special_tokens if g not in old_vocab]
        gene_ids_pretrained = np.array(old_vocab(pretrained_genes), dtype=int)
        # https://discuss.pytorch.org/t/expand-an-existing-embedding-and-linear-layer-nan-loss-value/55670/2
        # Retrieve pretrained weights
        vocab = Vocab(VocabPybind(pretrained_genes + new_genes, None))
        vocab.set_default_index(vocab["<pad>"])
        gene_ids = np.array(vocab(genes), dtype=int)

    tokenized_train = tokenize_and_pad_batch(
        train_data,
        gene_ids,
        max_len=config.max_seq_len,
        vocab=vocab,
        pad_token=config.pad_token,
        pad_value=config.pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=config.include_zero_gene,
        mod_type=mod_type if config.use_mod else None,
        vocab_mod=vocab_mod if config.use_mod else None,
    )
    tokenized_valid = tokenize_and_pad_batch(
        valid_data,
        gene_ids,
        max_len=config.max_seq_len,
        vocab=vocab,
        pad_token=config.pad_token,
        pad_value=config.pad_value,
        append_cls=True,
        include_zero_gene=config.include_zero_gene,
        mod_type=mod_type if config.use_mod else None,
        vocab_mod=vocab_mod if config.use_mod else None,
    )
    logger.info(
        f"train set number of samples: {tokenized_train['genes'].shape[0]}, "
        f"\n\t feature length: {tokenized_train['genes'].shape[1]}"
    )
    logger.info(
        f"valid set number of samples: {tokenized_valid['genes'].shape[0]}, "
        f"\n\t feature length: {tokenized_valid['genes'].shape[1]}"
    )

    device_number = 0
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{device_number}")
    else:
        device = torch.device("cpu")
    model_dict = torch.load(model_file)
    ntokens = len(vocab)  # size of vocabulary
    model = MultiOmicTransformerModel(
        ntokens,
        embsize,
        nhead,
        d_hid, 
        nlayers,
        vocab=vocab,
        dropout=config.dropout,
        pad_token=config.pad_token,
        pad_value=config.pad_value,
        do_mvc=config.GEPC,
        do_dab=config.DAR,
        use_batch_labels=config.use_batch_labels,
        num_batch_labels=num_batch_types,
        domain_spec_batchnorm=config.DSBN,
        n_input_bins=config.n_bins,
        ecs_threshold=config.ecs_thres,
        explicit_zero_prob=config.explicit_zero_prob,
        use_fast_transformer=config.fast_transformer,
        pre_norm=config.pre_norm,
        use_mod=config.use_mod,
        ntokens_mod=ntokens_mod if config.use_mod else None,
        vocab_mod=vocab_mod if config.use_mod else None,
    )

    with torch.no_grad():
        pretrained_emb_weights = model_dict['encoder.embedding.weight'][gene_ids_pretrained, :]
        model.encoder.embedding.weight.data[:len(pretrained_genes), :] = pretrained_emb_weights
        model.encoder.enc_norm.weight.data = model_dict['encoder.enc_norm.weight']

        if config.load_layers > 0:
            keep_keys = ['value_encoder'] + [f'transformer_encoder.layers.{i}' for i in range(config.load_layers)]
        else:
            keep_keys = []

        if args.pretrain:
            model_state_dict = model.state_dict()
            pretrained_state_dict = {}
            for key, value in model_dict.items():
                if key in model_state_dict:
                    if key in keep_keys:
                        pretrained_state_dict[key] = value

            model_state_dict.update(pretrained_state_dict)
            model.load_state_dict(model_state_dict)

    ntokens = len(vocab)

    model.to(device)
    print(model)

    if config.GEP and config.GEPC:
        criterion_gep_gepc = masked_mse_loss
    if config.CLS:
        criterion_cls = nn.CrossEntropyLoss()
    if config.DAR:
        criterion_dab = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(
        model.parameters(), lr=config.lr, eps=1e-4 if config.amp else 1e-8
    )
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1, gamma=config.schedule_ratio)
    scaler = torch.cuda.amp.GradScaler(enabled=config.amp)

    best_val_loss = float("inf")
    best_avg_bio = 0.0
    best_model = None
    define_wandb_metrcis()

    print("[5/6] Finetuning model...")
    for epoch in range(1, config.epochs + 1):
        epoch_start_time = time.time()
        train_data_pt, valid_data_pt = prepare_data(
            tokenized_train=tokenized_train, 
            tokenized_valid=tokenized_valid, 
            train_batch_labels=train_batch_labels,
            valid_batch_labels=valid_batch_labels,
            config=config,
            epoch=epoch,
            sort_seq_batch=config.per_seq_batch_sample)
        
        train_loader = prepare_dataloader(
            train_data_pt,
            batch_size=config.batch_size,
            shuffle=True,
            intra_domain_shuffle=False,
            drop_last=False,
            per_seq_batch_sample=config.per_seq_batch_sample
        )
        valid_loader = prepare_dataloader(
            valid_data_pt,
            batch_size=config.batch_size,
            shuffle=False,
            intra_domain_shuffle=False,
            drop_last=False,
            per_seq_batch_sample=config.per_seq_batch_sample
        )

        if config.do_train:
            train(
                model=model,
                loader=train_loader,
                vocab=vocab,
                criterion_gep_gepc=criterion_gep_gepc if config.GEP and config.GEPC else None,
                criterion_dab=criterion_dab if config.DAR else None,
                criterion_cls=criterion_cls if config.CLS else None,
                scaler=scaler,
                optimizer=optimizer,
                scheduler=scheduler,
                device=device,
                config=config,
                logger=logger,
                epoch=epoch,
            )
        val_loss = evaluate(
            model=model,
            loader=valid_loader,
            vocab=vocab,
            criterion_gep_gepc=criterion_gep_gepc if config.GEP and config.GEPC else None,
            criterion_dab=criterion_dab if config.DAR else None,
            criterion_cls=criterion_cls if config.CLS else None,
            device=device,
            config=config,
            epoch=epoch
        )
        elapsed = time.time() - epoch_start_time
        logger.info("-" * 89)
        logger.info(
            f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
            f"valid loss {val_loss:5.4f} | "
        )
        logger.info("-" * 89)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model = copy.deepcopy(model)
            best_model_epoch = epoch
            logger.info(f"Best model with score {best_val_loss:5.4f}")

        if epoch % config.save_eval_interval == 0 or epoch == config.epochs:
            logger.info(f"Saving model to {save_dir}")
            torch.save(best_model.state_dict(), save_dir / f"model_e{best_model_epoch}.pt")

            results,cell_embeddings = eval_testdata_save_embed(
                model = best_model,
                adata_t = adata_sorted if config.per_seq_batch_sample else adata,
                gene_ids = gene_ids,
                vocab = vocab,
                config = config,
                logger = logger,
                include_types=["cls"],
            )

            pd.DataFrame(cell_embeddings, index=adata.obs_names).to_csv(output_rna_dir, header=False)
            pd.DataFrame(cell_embeddings, index=adata.obs_names).to_csv(output_atac_dir, header=False)

            results["batch_umap"].savefig(
                save_dir / f"embeddings_batch_umap[cls]_e{best_model_epoch}.png", dpi=300
            )

            results["celltype_umap"].savefig(
                save_dir / f"embeddings_celltype_umap[cls]_e{best_model_epoch}.png", dpi=300
            )
            metrics_to_log = {"test/" + k: v for k, v in results.items()}
            metrics_to_log["test/batch_umap"] = wandb.Image(
                str(save_dir / f"embeddings_batch_umap[cls]_e{best_model_epoch}.png"),
                caption=f"celltype avg_bio epoch {best_model_epoch}",
            )

            metrics_to_log["test/celltype_umap"] = wandb.Image(
                str(save_dir / f"embeddings_celltype_umap[cls]_e{best_model_epoch}.png"),
                caption=f"celltype avg_bio epoch {best_model_epoch}",
            )
            metrics_to_log["test/best_model_epoch"] = best_model_epoch
            wandb.log(metrics_to_log)
            wandb.log({"avg_bio": results.get("avg_bio", 0.0)})

        scheduler.step()

    print("[6/6] Saving best model and results...")
    elapsed_time = time.time() - start_time
    with open(save_dir / 'run_info.yaml','w') as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": adata.shape[0]
        }, f)
    torch.save(best_model.state_dict(), save_dir / "best_model.pt")

if __name__ == "__main__":
    main(parse_args())