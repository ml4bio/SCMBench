#!/bin/bash

CUDA_VISIBLE_DEVICES=1 /mnt/nas/user/yixuan/miniconda3/envs/scMDC/bin/python /mnt/nas/user/yixuan/MMD_MA/run_MMD_MA.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/MMD_MA/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/MMD_MA/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/MMD_MA/run_info_gpu.json \