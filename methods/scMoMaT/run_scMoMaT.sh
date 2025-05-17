#!/bin/bash

CUDA_VISIBLE_DEVICES=1 /workspace/wangyixuan/.conda/envs/zgy_py3.8/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/scMoMaT/run_scMoMaT.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-feature /mnt/nas/user/yixuan/zgy/comp_effic/scMoMaT/feature_latent.csv \
    --train-dir /mnt/nas/user/yixuan/zgy/comp_effic/scMoMaT \
    --output-marker-feature /mnt/nas/user/yixuan/zgy/comp_effic/scMoMaT/marker_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/scMoMaT/run_info_gpu.json \
