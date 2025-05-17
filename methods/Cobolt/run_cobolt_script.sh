#!/bin/bash

dataset_name=10x-Multiome-Pbmc10k-small
for lr in 0.001 0.01 0.005; do
/ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_Cobolt.py \
    --input-joint-rna /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-RNA.h5ad \
    --input-joint-atac  /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-ATAC.h5ad \
    --output-joint /ailab/user/liuxinyuan/projects/scmbench/evaluation/test_cobolt/$dataset_name-$lr/$dataset_name-joint.csv \
    --lr $lr \
    -r /ailab/user/liuxinyuan/projects/scmbench/evaluation/test_cobolt/$dataset_name-$lr/$dataset_name-run-info.yaml
done