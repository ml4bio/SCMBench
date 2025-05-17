#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_py3.8/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/UnionCom/run_UnionCom.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/UnionCom/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/UnionCom/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/UnionCom/run_info_gpu.json \
    --device cuda