#!/bin/bash
#SBATCH --gpus=1
module load anaconda/2022.10
source activate scgpt

dataset_name=10x-Multiome-Pbmc10k-small
for lr in 0.001 0.01 0.005; do
/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scJoint.py \
    --input-rna /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-RNA.h5ad \
    --input-atac /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-ACTIVE.h5ad \
    --output-rna /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/$dataset_name-$lr/$dataset_name-rna.csv \
    --output-atac /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/$dataset_name-$lr/$dataset_name-atac.csv \
    --lr-stage1 $lr \
    --lr-stage3 $lr \
    -r /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/$dataset_name-$lr/$dataset_name-run-info.yaml
done