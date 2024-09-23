method=scGPT
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Integration_accuracy real: '$method'         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

for data_name in 10x-Multiome-Pbmc10k
do
  echo "Evaluation with $data_name"
  python integration_accuracy.py \
      -d ../../data/download/$data_name/$data_name-small-RNA_uni.h5ad ../../data/download/$data_name/$data_name-small-ATAC_uni.h5ad \
      -l ../../results/$method-output/$data_name/$data_name-small-rna.csv ../../results/$method-output/$data_name/$data_name-small-atac.csv \
      -o ../../results/$method-output/$data_name/$data_name-cell_integration_info.yaml 
  echo "=================="
done
