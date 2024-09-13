#!/bin/bash
#SBATCH --gpus=1
module load anaconda/2022.10
source activate scgpt

for data_name in 10x-Multiome-Pbmc10k-small
do
  echo "Evaluation with $data_name"
  /ailab/group/groups/aim/liuxinyuan/.conda/envs/scgpt/bin/python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/evaluation/scripts/cell_integration.py \
      -d /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name-RNA.h5ad /ailab/user/liuxinyuan/projects/scmbench/datasets/$data_name-ATAC.h5ad \
      -l /ailab/user/liuxinyuan/projects/scmbench/save/scGPT_scVI_30/$data_name-rna.csv /ailab/user/liuxinyuan/projects/scmbench/save/scGPT_scVI_30/$data_name-atac.csv \
      -o /ailab/user/liuxinyuan/projects/scmbench/save/scGPT_scVI_30/$data_name-cell_integration_info.yaml 
  echo "=================="
done
