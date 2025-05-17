#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_r/bin/Rscript /mnt/nas/user/yixuan/zgy/comp_effic/Seurat4/run_Seurat4_with_activity.R \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --activity /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ACTIVE.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/Seurat4/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/Seurat4/atac_latent.csv \
    --output-activity /mnt/nas/user/yixuan/zgy/comp_effic/Seurat4/activity_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/Seurat4/run_info_cpu.json \
