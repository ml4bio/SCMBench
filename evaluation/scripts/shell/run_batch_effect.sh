# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Integration: scMDC         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scMDC #scVI iNMF scMDC

# # for data_name in 10x-Multiome-Pbmc10k Chen-2019 Ma-2020 Yao-2021
# for data_name in nips_10x
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
# echo '+++++         Batch_effect: scVI_simu        +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI

# # for num_cell in 1000 5000 10000
# for num_cell in 1000
# do
#   # for num_batch in 2 3 5
#   for num_batch in 2
#   do
#     # for num_gene in 1000 3000
#     for num_gene in 1000
#     do
#       data_name='simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene
#       echo "Running exp with $data_name"
#       CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/batch_effect_metrics.py \
#             -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-RNA.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-ATAC.h5ad \
#             -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-RNA.csv   /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-ATAC.csv  \
#             -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated/$data_name-batch_effect_info.yaml 
#       echo "=================="
#     done
#   done
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
#         CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/batch_effect_metrics.py \
#           -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-RNA.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-ATAC.h5ad \
#           -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-rna.csv   /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-atac.csv  \
#           -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-batch_effect_info.yaml 
#         echo "=================="
#       done
#     done
#   done
# done


method=scFoundation #scVI iNMF scMDC UCE Geneformer scFoundation scGPT-zero scGPT

echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Integration_batch effect real: '$method'         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# m=kidney #bc
# for data_name in 10x-Multiome-Pbmc10k Chen-2019 Ma-2020 Yao-2021 Muto-2021 Ma-2020-small
# for data_name in Ma-2020-sampled-small Ma-2020-small 
for data_name in Muto-2021-sampled-small Muto-2021-small
do
  echo "Running exp with $data_name"
  CUDA_VISIBLE_DEVICES=0 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/batch_effect_metrics.py \
        -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA_uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC_uni.h5ad\
        -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-atac.csv \
        -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-batch_effect_info.yaml 
  echo "=================="
done