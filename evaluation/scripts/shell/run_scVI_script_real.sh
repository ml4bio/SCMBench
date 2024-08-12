# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: scVI         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI

# # for data_name in 10x-Multiome-Pbmc10k Chen-2019 Muto-2021 Ma-2020 Yao-2021
# for data_name in nips_challenge_cite_paired 
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#         -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
#   echo "=================="
# done


# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: scVI_nips         +++++'
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
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/$data_name-${data_type,,}.csv  \
#         -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/run_info.yaml
#   echo "=================="
#   data_type='RNA'
#   echo "Running exp with $dir_name $data_type"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name-$data_type.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/$data_name-${data_type,,}.csv  \
#         -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/run_info.yaml
#   echo "=================="
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Multiomics: scVI_nips_threeomics         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI

# for data_name in rna adt atac
# do
#   dir_name='nips_challenge_triple_data_unpaired'
#   echo "Running exp with $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#         --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$dir_name/$data_name.h5ad \
#         --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/$data_name.csv  \
#         -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$dir_name/run_info.yaml
#   echo "=================="

# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         simulate_batch_effect: scVI         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI
# # for num_cell in 1000 5000 10000
# for num_cell in 2000 5000 10000
# do
#   for num_batch in 5
#   do
#     for num_gene in 5000
#     do
#       for effect in 20 50 80
#       do
#         for data_type in RNA ATAC
#         do
#           data_name='0910_grn_1139_simulated_num_cell_'$num_cell'_num_batch_'$num_batch'_num_gene_'$num_gene'-noised-effect-'$effect'-'$data_type
#           echo "Running exp with $data_name"
#           CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#                 --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/simulated_var_cell_effect/$data_name.h5ad \
#                 --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name.csv  \
#                 -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/preprocessed_simulated_noised_effect/$data_name-run_info.yaml
#           echo "=================="
#         done
#       done
#     done
#   done
# done

# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         10X_nips: scVI         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=scVI
# data_name=nips_10x

# echo "Running exp with $data_name -ATAC"
# CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#       --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-ATAC.h5ad \
#       --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#       -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
# echo "=================="

# echo "Running exp with $data_name -RNA"
# CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
#       --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/seurat_10x/$data_name/$data_name-RNA.h5ad \
#       --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#       -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
# echo "=================="

echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         small: scVI         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

method=scVI
# for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# for data_name in Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
for data_name in Muto-2021-sampled-small
do
      echo "Running exp with $data_name -ATAC"
      CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
            --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad \
            --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
            -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
      echo "=================="

      echo "Running exp with $data_name -RNA"
      CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/run_scVI.py \
            --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
            --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
            -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
      echo "=================="
done