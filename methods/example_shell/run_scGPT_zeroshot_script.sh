#!/bin/bash

dataset_name=10x-Multiome-Pbmc10k-small
python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scGPT_zeroshot.py \
--input-path /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-RNA.h5ad \
--model-path /ailab/user/liuxinyuan/projects/foundation_models/scGPT/blood \
--output-path /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/scgpt_zero_rna.csv

python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/run/run_scGPT_zeroshot.py \
--input-path /ailab/user/liuxinyuan/projects/scmbench/datasets/$dataset_name/$dataset_name-ACTIVE.h5ad \
--model-path /ailab/user/liuxinyuan/projects/foundation_models/scGPT/blood \
--output-path /ailab/user/liuxinyuan/projects/scmbench/evaluation/test/scgpt_zero_atac.csv
