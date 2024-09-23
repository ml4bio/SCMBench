#!/bin/bash
#SBATCH --gpus=1
module load anaconda/2022.10
source activate perturb

data_name=10x-Multiome-Pbmc10k-small
python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_FM_scVI.py \
--input-rna /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
--input-atac /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
--rna-pre /ailab/user/liuxinyuan/projects/scmbench/evaluation/scGPT-zero-output/$data_name/$data_name-human-rna.csv \
--atac-pre /ailab/user/liuxinyuan/projects/scmbench/evaluation/scGPT-zero-output/$data_name/$data_name-human-atac.csv \
--output-path /ailab/user/liuxinyuan/projects/scmbench/evaluation/test_scgpt_scvi/


