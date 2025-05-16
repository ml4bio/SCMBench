echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Multiomics: Deepmaps         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

method=Deepmaps

# for data_name in Chen-2019 Ma-2020 Yao-2021 10x-Multiome-Pbmc10k
# for data_name in 10x-Multiome-Pbmc10k
# for data_name in 10x-Multiome-Pbmc10k-small
# for data_name in Ma-2020-small-10000
# for data_name in Ma-2020-small-4000
# for data_name in 10x-Multiome-Pbmc10k-small 
# do
#   echo "Running exp with Human $data_name"
#   Rscript /mnt/nas/user/yixuan/deepmaps-master/run_Deepmaps.R \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
#         --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#         --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#         --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml \
#         --genome GRCh38 
#   echo "=================="
# done

# for data_name in Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small

for data_name in Ma-2020-small
do
  echo "Running exp with Mouse $data_name"
  Rscript /mnt/nas/user/yixuan/deepmaps-master/run_Deepmaps.R \
        --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
        --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
        --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
        --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
        --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml \
        --genome GRCm38 
  echo "=================="
done


# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: Deepmaps_simu        +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=Deepmaps

# for num_cell in 1000 5000 10000
# do
#   for num_batch in 2 3 5
#   do
#     for num_gene in 1000 3000
#     do
#       data_name='simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene
#       echo "Running exp with $data_name"
#       Rscript /mnt/nas/user/yixuan/deepmaps-master/run_Deepmaps.R \
#             --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-RNA.h5ad \
#             --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-ATAC.h5ad \
#             --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-rna.csv  \
#             --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-atac.csv  \
#             --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-run_info.yaml
#       echo "=================="
#     done
#   done
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         simulate_batch_effect: Deepmaps         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=Deepmaps
# # for num_cell in 1000 5000 10000
# for num_cell in 2000 5000 10000
# do
#   for num_batch in 5
#   do
#     for num_gene in 5000
#     do
#       for effect in 20 50 80
#       do
#         data_name='0910_grn_1139_simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene'-noised-effect-'$effect
#         echo "Running exp with $data_name"
#         Rscript /mnt/nas/user/yixuan/deepmaps-master/run_Deepmaps.R \
#               --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-RNA.h5ad \
#               --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-ATAC.h5ad \
#               --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-rna.csv  \
#               --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-atac.csv  \
#               --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-run_info.yaml 
#         echo "=================="
#       done
#     done
#   done
# done

# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         10X_nips: simulate_batch_effect         +++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=Deepmaps
# data_name=nips_10x

# echo "Running exp with $data_name"

#   Rscript /mnt/nas/user/yixuan/deepmaps-master/run_Deepmaps.R \
#       --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-RNA_dense.h5ad \
#       --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-ATAC_dense.h5ad \
#       --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#       --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#       --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-run_info.yaml 
# echo "=================="
