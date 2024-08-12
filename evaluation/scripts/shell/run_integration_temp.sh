echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Integration_batch effect real: scMDC         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

method=Deepmaps #scVI iNMF scMDC

# for data_name in 10x-Multiome-Pbmc10k Chen-2019 Yao-2021
# for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small 
for data_name in Ma-2020-small
do
  echo "Running exp with $data_name"
  CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
        -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA_uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC_uni.h5ad\
        -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv \
        -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/cell_integration_info.yaml \
        -p
  echo "=================="
done
# for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
# # for data_name in Muto-2021-small Yao-2021-small
# do
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
#         -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA_uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC_uni.h5ad\
#         -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv \
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/cell_integration_info.yaml 
#   echo "=================="
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Integration: scVI_nips         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI

# for data_name in cite multiome
# do
#   if [ $data_name == 'cite' ] 
#   then  
#     data_type='ADT'
#   else
#     data_type='ATAC'
#   fi  
#   dir_name='nips_challenge_'$data_name'_paired'
#   echo "Running exp with $dir_name $data_type"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
#         -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type-uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-RNA-uni.h5ad \
#         -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/$data_name-${data_type,,}.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/$data_name-rna.csv  \
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/cell_integration_info.yaml 
#   echo "=================="
# done


# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: scVI_simu        +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI

# for num_cell in 1000 5000 10000
# do
#   for num_batch in 2 3 5
#   do
#     for num_gene in 1000 3000
#     do
#       data_name='simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene
#       echo "Running exp with $data_name"
#       CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
#             -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-RNA-uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-ATAC-uni.h5ad \
#             -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-RNA.csv   /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-ATAC.csv  \
#             -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-cell_integration_info.yaml 
#       echo "=================="
#     done
#   done
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Integration_batch effect real: scMDC         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=iNMF #scVI iNMF scMDC

# # for data_name in 10x-Multiome-Pbmc10k Chen-2019 Ma-2020 Yao-2021
# for data_name in nips_10x
# do
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
#         -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-RNA_uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-ATAC_uni.h5ad\
#         -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv \
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/cell_integration_info.yaml 
#   echo "=================="
# done


# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         simulate_batch_effect: Batch_effect         +++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=iNMF

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
#         CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
#           -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-RNA-uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-ATAC-uni.h5ad \
#           -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-rna.csv   /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-atac.csv  \
#           -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-cell_integration_info.yaml 
#         echo "=================="
#       done
#     done
#   done
# done
