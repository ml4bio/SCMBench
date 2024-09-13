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
  CUDA_VISIBLE_DEVICES=0 python /ailab/user/liuxinyuan/projects/scmbench/SCMBench/evaluation/scripts/batch_effect_metrics.py \
        -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA_uni.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC_uni.h5ad\
        -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-atac.csv \
        -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-batch_effect_info.yaml 
  echo "=================="
done