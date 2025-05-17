#!/usr/bin/env bash

CUDA_VISIBLE_DEVICES=0 /workspace/wangyixuan/.conda/envs/zgy_py3.8/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/MOFA/run_MOFA_gpu.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/MOFA/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/MOFA/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/MOFA/run_info_gpu.json \
