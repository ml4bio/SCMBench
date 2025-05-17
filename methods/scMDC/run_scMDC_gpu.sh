#!/bin/bash

CUDA_VISIBLE_DEVICES=1 /workspace/wangyixuan/.conda/envs/zgy_py3.8/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/scMDC/run_scMDC.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/scMDC/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/scMDC/atac_latent.csv \
    --save_dir /mnt/nas/user/yixuan/zgy/comp_effic/scMDC \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/scMDC/run_info_gpu.json \
    --device cuda \
    --ae_weight_file /mnt/nas/user/yixuan/zgy/comp_effic/scMDC/AE_weights_1.pth.tar