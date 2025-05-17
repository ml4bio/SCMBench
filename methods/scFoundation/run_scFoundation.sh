#!/bin/bash

/mnt/nas/user/yixuan/miniconda3/envs/zgy_py3.10/bin/python /mnt/nas/user/yixuan/scFoundation/model/run_scFoundation.py \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad \
    --output-rna /mnt/nas/user/yixuan/scFoundation/model/rna_latent.csv \
    --run-info /mnt/nas/user/yixuan/scFoundation/model/run_info_gpu.json \
    --gene_list /mnt/nas/user/yixuan/scFoundation/model/OS_scRNA_gene_index.19264.tsv \
    --ckpt_path /mnt/nas/user/yixuan/scFoundation/model/models/models.ckpt \
    --device cuda