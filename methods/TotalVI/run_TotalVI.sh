#!/bin/bash

/workspace/wangyixuan/.conda/envs/scVI/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/TotalVI/run_TotalVI.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-protein /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --train-dir /mnt/nas/user/yixuan/zgy/comp_effic/TotalVI \
    --output-feature /mnt/nas/user/yixuan/zgy/comp_effic/TotalVI/feature_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/TotalVI/run_info_gpu.json \
