#!/bin/bash

/workspace/wangyixuan/.conda/envs/scgpt/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/Geneformer/run_Geneformer.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/Geneformer/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/Geneformer/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/Geneformer/run_info_gpu.json \