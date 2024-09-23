# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: iNMF         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=iNMF
# unset LD_LIBRARY_PATH
# # for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# # for data_name in 10x-Multiome-Pbmc10k Chen-2019 Muto-2021 Ma-2020 Yao-2021 
# # for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
# for data_name in Muto-2021-sampled-small
# do
#   echo "Running exp with $data_name"
#   Rscript /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_iNMF.R \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
#         --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#         --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#         --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml \

#   echo "=================="
# done



# for data_name in Muto-2021-sampled-small
# do
#   echo "Running exp with $data_name"
#   Rscript /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_iNMF.R \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
#         --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#         --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#         --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml \

#   echo "=================="
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: iNMF_simu        +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=iNMF
# unset LD_LIBRARY_PATH
# for num_cell in 1000 5000 10000
# do
#   for num_batch in 2 3 5
#   do
#     for num_gene in 1000 3000
#     do
#       data_name='simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene
#       echo "Running exp with $data_name"
#       Rscript /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_iNMF.R \
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
# echo '+++++         simulate_batch_effect: iNMF         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=iNMF
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
#         Rscript /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_iNMF.R \
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

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         simulate_batch_effect: iNMF         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

method=iNMF
# for num_cell in 1000 5000 10000
for num_cell in 5000
do
  for num_batch in 3
  do
    for num_gene in 3000
    do
      for effect in 0 0.5 1 1.5 2 3 5 10
      do
        data_name='simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene'_effect_'$effect
        echo "Running exp with $data_name"
        Rscript /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_iNMF.R \
              --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulate_new_for_test/$data_name-RNA.h5ad \
              --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulate_new_for_test/$data_name-ATAC.h5ad \
              --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/simulate_new_for_test/$data_name-rna.csv  \
              --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/simulate_new_for_test/$data_name-atac.csv  \
              --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/simulate_new_for_test/$data_name-run_info.yaml 
        echo "=================="
      done
    done
  done
done

# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         10X_nips: simulate_batch_effect         +++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=iNMF
# data_name=nips_10x

# echo "Running exp with $data_name"

#   Rscript /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_iNMF.R \
#       --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-RNA_dense.h5ad \
#       --input-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-ATAC_dense.h5ad \
#       --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#       --output-atac /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#       --run-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-run_info.yaml 
# echo "=================="
