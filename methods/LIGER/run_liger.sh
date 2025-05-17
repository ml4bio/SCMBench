#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_r/bin/Rscript /mnt/nas/user/yixuan/zgy/comp_effic/liger/run_LIGER.R \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ACTIVE.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/liger/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/liger/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/liger/run_info_cpu.json \
