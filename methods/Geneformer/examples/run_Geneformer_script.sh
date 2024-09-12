
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         small: Geneformer         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=Geneformer
# # for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# # for data_name in Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# # for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
# for data_name in 10x-Multiome-Pbmc10k-small
# do

#       echo "Running exp with $data_name -RNA"
#       CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Geneformer/examples/run_Geneformer.py \
#             --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
#             --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  \
#             -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
#       echo "=================="
# done



# echo '++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         small: Geneformer         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++'

# method=Geneformer
# # for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# # for data_name in Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# # for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
# for data_name in 10x-Multiome-Pbmc10k-small
# do
#       echo "Running exp with $data_name -ATAC Activity"
#       CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Geneformer/examples/run_Geneformer_ACTIVE.py \
#             --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ACTIVE.h5ad \
#             --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv  \
#             -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
#       echo "=================="
# done

echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         small: Geneformer         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

method=Geneformer
# for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# for data_name in Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small Muto-2021-small Yao-2021-small
# for data_name in Muto-2021-sampled-small Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
for data_name in 10x-Multiome-Pbmc10k-small
do
      echo "Running exp with $data_name -ATAC Activity"
      CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Geneformer/examples/run_Geneformer_ACTIVE.py \
            --input-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ACTIVE-con.h5ad \
            --output-rna /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-con.csv  \
            -r /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/run_info.yaml
      echo "=================="
done