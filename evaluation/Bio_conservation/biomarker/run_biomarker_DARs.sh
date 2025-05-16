echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Biomarker                    +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

method = scVI GLUE LIGER TotalVI UnionCom scMoMaT MMD_MA Deepmaps Cobolt scMDC bindsc Pamona PCA scJoint seurat4 seurat5 iNMF MOFA Harmony UCE FM_scvi
for data_name in Ma-2020-batch-53-small
do
  echo "Running exp with $data_name"
  python biomarker.py \
        -d $data_name/$data_name-ATAC.h5ad \
        -l ../../../Saved_embeddings \
        -m $method \
        --mode atac \
        -o ../../eval_metrics/biomarker/DARs/$data_name
  echo "=================="
done