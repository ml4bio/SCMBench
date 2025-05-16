#!/bin/bash
#SBATCH --gpus=1
module load anaconda/2022.10
source activate scgpt

/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_uce.py \
    --adata_path /ailab/user/liuxinyuan/projects/scmbench/datasets/Chen-2019-small/Chen-2019-small-RNA.h5ad \
    --output_dir /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/UCE-zero-output/Chen-2019-small/ \
    --species mouse