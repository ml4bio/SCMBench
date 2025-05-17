#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_r/bin/Rscript /mnt/nas/user/yixuan/zgy/comp_effic/Seurat5/run_Seurat5.R \
    --input-bridge-rna /mnt/nas/user/yixuan/SCMBench/data/download/Muto-2021-batch-1-small/Muto-2021-batch-1-small-RNA.h5ad \
    --input-bridge-atac /mnt/nas/user/yixuan/SCMBench/data/download/Muto-2021-batch-1-small/Muto-2021-batch-1-small-ATAC.h5ad \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/Muto-2021-batch-1-small/Muto-2021-batch-1-small-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/Muto-2021-batch-1-small/Muto-2021-batch-1-small-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/Seurat5/rna_latent_Muto.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/Seurat5/atac_latent_Muto.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/Seurat5/run_info_Muto.json \
    --device cpu
