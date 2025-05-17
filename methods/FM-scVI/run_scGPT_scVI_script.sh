

data_name=10x-Multiome-Pbmc10k-small
python /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/SCMBench/methods/run_gpt_scvi.py \
--input-rna /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-RNA.h5ad \
--input-atac /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/datasets/$data_name/$data_name-ATAC.h5ad \
--rna-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scGPT-zero-output/$data_name/$data_name-human-rna.csv \
--atac-pre /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scGPT-zero-output/$data_name/$data_name-human-atac.csv \
--output-path /ailab/user/chenpengan/liuxinyuan/new_liuxinyuan/scmbench/evaluation/scgpt_scvi/$data_name/
