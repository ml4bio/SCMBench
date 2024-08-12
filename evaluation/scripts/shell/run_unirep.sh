# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         UNI: ATAC         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# for data_name in 10x-Multiome-Pbmc10k Chen-2019 Muto-2021 Ma-2020 Yao-2021
# do
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/atac_unirep.py \
#         -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad\
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC_uni.h5ad
#   echo "=================="
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         UNI: RNA         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# for data_name in 10x-Multiome-Pbmc10k Chen-2019 Muto-2021 Ma-2020 Yao-2021
# do
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rna_unirep.py \
#         -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad\
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA_uni.h5ad
#   echo "=================="
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         UNI: nips         +++++'
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
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/atac_unirep.py \
#         -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type.h5ad\
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type-uni.h5ad
#   echo "=================="
#   data_type='RNA'
#   echo "Running exp with $dir_name $data_type"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rna_unirep.py \
#         -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type.h5ad\
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type-uni.h5ad
#   echo "=================="
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         UNI: simu        +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scMDC

# for num_cell in 1000 5000 10000
# do
#   for num_batch in 2 3 5
#   do
#     for num_gene in 1000 3000
#     do
#       data_name='simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene
#       echo "Running exp with $data_name"
#       CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/atac_unirep.py \
#             -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-ATAC.h5ad \
#             -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-ATAC-uni.h5ad 
#       echo "=================="
#       CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rna_unirep.py \
#             -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-RNA.h5ad \
#             -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/preprocessed_simulated/$data_name-RNA-uni.h5ad
#       echo "=================="
#     done
#   done
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: iNMF         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scMDC

# # for data_name in 10x-Multiome-Pbmc10k Chen-2019 Ma-2020 Yao-2021
# for data_name in Muto-2021

# do
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py \
#         -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad\
#         -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv \
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/cell_integration_info.yaml
#   echo "=================="
# done


#Deepmaps

# f1=/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-RNA.h5ad
# f11=/mnt/nas/user/yixuan/deepmaps-master/debug-output/rna.csv
# f2=/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/10x-Multiome-Pbmc10k/10x-Multiome-Pbmc10k-ATAC.h5ad
# f21=/mnt/nas/user/yixuan/deepmaps-master/debug-output/atac.csv
# CUDA_VISIBLE_DEVICES=1 python -u /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/cell_integration.py -d $f1 $f2 -l $f11 $f21 -o /mnt/nas/user/yixuan/deepmaps-master/debug-output/cell_integration_info.yaml


# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         10X_nips: scVI         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI
# data_name=nips_10x

# echo "Running exp with $data_name -ATAC"
# CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/atac_unirep.py \
#       -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-ATAC.h5ad \
#       -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-ATAC_uni.h5ad 
# echo "=================="

# echo "Running exp with $data_name -RNA"
# CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rna_unirep.py \
#       -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-RNA.h5ad\
#       -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-RNA_uni.h5ad
# echo "=================="

echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '++++++++++++         unirep         ++++++++++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
for data_name in Muto-2021-sampled-small
do
      echo "Running exp with $data_name -ATAC"
      CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/atac_unirep.py \
            -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
            -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC_uni.h5ad 
      echo "=================="

      echo "Running exp with $data_name -RNA"
      CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rna_unirep.py \
            -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad\
            -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA_uni.h5ad
      echo "==================="
done
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         simulate_batch_effect: Batch_effect         +++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

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
#         CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/atac_unirep.py \
#             -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-ATAC.h5ad \
#             -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-ATAC-uni.h5ad 
#         echo "=================="
#         CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rna_unirep.py \
#             -i /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-RNA.h5ad \
#             -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name-RNA-uni.h5ad
#         echo "=================="
#       done
#     done
#   done
# done
