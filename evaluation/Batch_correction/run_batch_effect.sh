method=scFoundation #scVI iNMF scMDC UCE Geneformer scFoundation scGPT-zero scGPT

echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Integration_batch effect real: '$method'         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


for data_name in Muto-2021-small
do
  echo "Running exp with $data_name"
  CUDA_VISIBLE_DEVICES=0 python batch_effect_metrics.py \
        -d $data_name-RNA_uni.h5ad $data_name-ATAC_uni.h5ad \
        -l $data_name$m-rna.csv  $data_name$m-atac.csv \
        -o $data_name$m-batch_effect_info.yaml 
  echo "=================="
done