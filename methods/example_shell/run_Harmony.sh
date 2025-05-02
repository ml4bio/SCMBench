#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1 # change this the number of GPUs needed
#SBATCH --cpus-per-task=16 # don't change this
#SBATCH --output log_harmony.txt
#SBATCH -p optimal
#SBATCH -A optimal
# nvidia-smi

module load anaconda/2022.10
source activate /ailab/user/chenpengan/.conda/envs/zgy_py3.8

data_name='10x-Multiome-Pbmc10k-small'
echo $data_name
python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/run_Harmony.py \
  --input-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
  --input-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ACTIVE.h5ad \
  --output-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/harmony/10x-Multiome-Pbmc10k/$data_name-RNA_latent2.csv  \
  --output-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/harmony/10x-Multiome-Pbmc10k/$data_name-ACTIVE_latent2.csv  \
  -r /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/harmony/10x-Multiome-Pbmc10k/$data_name-run_info2.yaml

