#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_r/bin/Rscript /mnt/nas/user/yixuan/zgy/comp_effic/iNMF/run_iNMF.R \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ACTIVE.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/iNMF/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/iNMF/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/iNMF/run_info_cpu.json \
    --device cpu
