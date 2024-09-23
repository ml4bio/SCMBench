#!/bin/bash
#SBATCH --gpus=1
module load anaconda/2022.10
source activate scgpt

data_name=10x-Multiome-Pbmc10k-small

/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scGPT_ft.py \
--input-rna /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
--input-atac /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
--model-path /ailab/user/liuxinyuan/projects/foundation_models/scGPT/whole_human \
--output-dir /ailab/user/liuxinyuan/projects/scmbench/evaluation/test_scgpt \
--pretrain False

/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scGPT_ft.py \
--input-rna /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
--input-atac /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
--model-path /ailab/user/liuxinyuan/projects/foundation_models/scGPT/whole_human \
--output-dir /ailab/user/liuxinyuan/projects/scmbench/evaluation/test_scgpt \
--pretrain True
