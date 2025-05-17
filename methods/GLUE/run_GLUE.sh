#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_py3.8/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/run_GLUE.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --prior /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/guidance.graphml.gz \
    --train-dir /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/ \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/atac_latent.csv \
    --output-feature /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/feature_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/GLUE/run_info_gpu.json
