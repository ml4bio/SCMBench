#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1 # change this the number of GPUs needed
#SBATCH --cpus-per-task=16 # don't change this
#SBATCH --output /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/log.txt
#SBATCH -p optimal
#SBATCH -A optimal 
# nvidia-smi

module load anaconda/2022.10
source activate gpt

# data_name=10x-Multiome-Pbmc10k-small
# python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/methods/run_gpt_scvi.py \
# --input-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
# --input-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
# --rna-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scGPT-zero-output/$data_name/$data_name-human-rna.csv \
# --atac-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scGPT-zero-output/$data_name/$data_name-human-atac.csv \
# --output-path /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/

# echo "=================="
# echo Evaluation with $data_name
# python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/tutorials/metrics/integration_accuracy.py \
#     -d /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA_uni.h5ad /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC_uni.h5ad \
#     -l /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/rna_embeddings.csv /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/atac_embeddings.csv \
#     -o /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/cell_integration_info.yaml 
# echo "=================="

# for data_name in Chen-2019-small Ma-2020-small; do
for data_name in TC; do
python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/methods/run_gpt_scvi.py \
--input-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
--input-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
--rna-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/UCE-zero-output/$data_name/rna_embeddings.csv \
--atac-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/UCE-zero-output/$data_name/atac_embeddings.csv \
--output-path /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/

echo "=================="
echo Evaluation with $data_name
python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/tutorials/metrics/integration_accuracy.py \
    -d /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA_uni.h5ad /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC_uni.h5ad \
    -l /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/rna_embeddings.csv /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/atac_embeddings.csv \
    -o /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/cell_integration_info.yaml 
echo "=================="

done

# for data_name in 10x-Multiome-Pbmc10k-small; do

# python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/methods/run_gpt_scvi.py \
# --input-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
# --input-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
# --rna-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scGPT-zero-output/$data_name/rna_embeddings.csv \
# --atac-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scGPT-zero-output/$data_name/atac_embeddings.csv \
# --output-path /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/

# echo "=================="
# echo Evaluation with $data_name
# python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/tutorials/metrics/integration_accuracy.py \
#     -d /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA_uni.h5ad /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC_uni.h5ad \
#     -l /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/rna_embeddings.csv /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/atac_embeddings.csv \
#     -o /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/cell_integration_info.yaml 
# echo "=================="

# done
