#!/bin/bash
#SBATCH --gpus=1
module load anaconda/2022.10
source activate scgpt

dataset_name=10x-Multiome-Pbmc10k-small
/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scGPT_zeroshot.py \
--input-path /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-RNA.h5ad \
--model-path /ailab/user/liuxinyuan/projects/foundation_models/scGPT/blood \
--output-path /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/scgpt_zero_rna.csv

/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scGPT_zeroshot.py \
--input-path /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-ACTIVE.h5ad \
--model-path /ailab/user/liuxinyuan/projects/foundation_models/scGPT/blood \
--output-path /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/scgpt_zero_atac.csv
