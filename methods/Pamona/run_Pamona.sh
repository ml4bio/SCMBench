#!/bin/bash

/workspace/wangyixuan/.conda/envs/zgy_py3.8/bin/python /mnt/nas/user/yixuan/zgy/comp_effic/Pamona/run_Pamona.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/Pamona/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/Pamona/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/Pamona/run_info_cpu.json
