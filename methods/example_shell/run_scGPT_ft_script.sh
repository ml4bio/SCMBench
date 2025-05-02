#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1 # change this the number of GPUs needed
#SBATCH --cpus-per-task=16 # don't change this
#SBATCH --output /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/log.txt
#SBATCH -p optimal
#SBATCH -A optimal 
# nvidia-smi

module load anaconda/2022.10
source activate gpt

data_name=10x-Multiome-Pbmc10k-small

for nlayers in 4; do
for pretrain in False; do
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Finetune FM: 'scgpt'         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/methods/run_scGPT_ft2.py \
--input-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
--input-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
--model-path /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/foundation_models/scGPT/whole_human \
--output-dir /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_ft_results_all/${data_name}_${nlayers}_pre${pretrain} \
--pretrain $pretrain
--load_layers $nlayers
--nlayers $nlayers

echo "=================="
echo Evaluation with $data_name
python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/tutorials/metrics/integration_accuracy.py \
    -d /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA_uni.h5ad /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC_uni.h5ad \
    -l /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_ft_results_all/${data_name}_${nlayers}_pre${pretrain}/rna_embeddings.csv /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_ft_results_all/${data_name}_${nlayers}_pre${pretrain}/atac_embeddings.csv \
    -o /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_ft_results_all/${data_name}_${nlayers}_pre${pretrain}/cell_integration_info.yaml 
echo "=================="

done;
done

