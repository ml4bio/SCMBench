import os
import sys
from pathlib import Path
import warnings

import scanpy as sc
import scib
import pandas as pd
import numpy as np
import sys

sys.path.insert(0, "../")

import scgpt as scg
import matplotlib.pyplot as plt
import anndata

plt.style.context('default')
warnings.simplefilter("ignore", ResourceWarning)

# model_dir = Path("../save/scGPT_bc")
model_dir = Path("../save/scGPT_human")

# dataset_dir_name='10x-Multiome-Pbmc10k-small'
# ls=['Muto-2021-batch-1-small','Muto-2021-batch-2-small','Muto-2021-batch-3-small','Muto-2021-batch-4-small','Muto-2021-batch-5-small','Muto-2021-sampled-small','Muto-2021-small']
# ls=['10x-Multiome-Pbmc10k-small']
ls=['Yao-2021-small']
for dataset_dir_name in ls:
    print('starts ATAC',dataset_dir_name)
    smaple_data_path='/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/'+dataset_dir_name+'/'+dataset_dir_name+'-ACTIVE.h5ad'

    # smaple_data_path = '../data/Kim2020_Lung.h5ad'
    adata = sc.read_h5ad(smaple_data_path)
    adata.var['gene_name']=list(adata.var.index)

    gene_col = "gene_name"
    cell_type_key = "cell_type"
    batch_key = "domain"

    celltype_id_labels = adata.obs[cell_type_key].astype("category").cat.codes.values
    adata = adata[celltype_id_labels >= 0]

    org_adata = adata.copy()

    embed_adata = scg.tasks.embed_data(
        adata,
        model_dir,
        gene_col=gene_col,
        batch_size=64,
    )
    # attach the cell embedding to the original adata

    emb=pd.DataFrame(embed_adata.obsm['X_scGPT'],index=adata.obs_names)
    emb.to_csv('/mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/scGPT-zero-output/'+dataset_dir_name+'/'+dataset_dir_name+'bc-atac.csv',header=False)