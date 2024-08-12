
# # m=kidney #bc kidney
# method=UCE  #scVI iNMF scMDC Geneformer scFoundation scGPT-zero scGPT UCE
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '+++++         Integration_trajectory real: '$method'         +++++'
# echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# # 10X-multiome (HSPC)
# # Chen. (OliI)
# # Ma-2020-batch-53-small (Basal)
# # Ma-2020-batch-54-small (Basal)
# # Ma-2020-batch-55-small (Basal)
# # Ma-2020-batch-56-small (Basal)
# # Ma-2020-sampled-small (Basal)
# # Ma-2020-small (Basal)
# # Muto-2021-batch-1-small (MES_FIB)
# # Muto-2021-batch-2-small (MES_FIB)
# # Muto-2021-batch-3-small (MES_FIB)
# # Muto-2021-batch-4-small (MES_FIB)
# # Muto-2021-batch-5-small (MES_FIB)
# # Muto-2021-sampled-small (MES_FIB)
# # Muto-2021-small (MES_FIB)
# # Yao-2021-small (NP)
# # paired_root=("HSPC" "OliI" "Basal" "Basal" "Basal" "Basal" "Basal" "Basal")

# pr=HSPC
# for data_name in 10x-Multiome-Pbmc10k-small
# # pr=OliI
# # for data_name in Chen-2019-small
# # pr=Basal
# # for data_name in Ma-2020-batch-53-small Ma-2020-batch-54-small Ma-2020-batch-55-small Ma-2020-batch-56-small Ma-2020-small Ma-2020-sampled-small
# # pr=MES_FIB
# # for data_name in Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
# # pr=NP
# # for data_name in Yao-2021-small

# # for data_name in Ma-2020-batch-55-small Chen-2019-small
# # for data_name in Ma-2020-batch-53-small Ma-2020-batch-54-small Ma-2020-batch-55-small Ma-2020-batch-56-small Ma-2020-small Ma-2020-sampled-small
# # for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small 
# # for data_name in Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small 
# do
#   echo "Running trajectory with $method $data_name"
#   CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/trajectory_conservation_metrics.py \
#         -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad \
#         -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-rna.csv \
#         -r $pr  \
#         -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/trajectory_info_$m.yaml
#   echo "=================="
# done




# m=kidney #bc kidney
# for method in bindsc Cobolt Deepmaps Geneformer GLUE iNMF LIGER MMD_MA MOFA Pamona PCA scFoundation scGPT scGPT-zero scJoint scMDC scMoMaT scVI seurat4 seurat5 TotalVI UCE UnionCom
for method in scMDC
do 
  # method=Deepmaps  #scVI iNMF scMDC Geneformer scFoundation scGPT-zero scGPT UCE
  echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  echo '+++++         Integration_trajectory real: '$method'         +++++'
  echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  # 10X-multiome (HSPC)
  # Chen. (OliI)
  # Ma-2020-batch-53-small (Basal)
  # Ma-2020-batch-54-small (Basal)
  # Ma-2020-batch-55-small (Basal)
  # Ma-2020-batch-56-small (Basal)
  # Ma-2020-sampled-small (Basal)
  # Ma-2020-small (Basal)
  # Muto-2021-batch-1-small (MES_FIB)
  # Muto-2021-batch-2-small (MES_FIB)
  # Muto-2021-batch-3-small (MES_FIB)
  # Muto-2021-batch-4-small (MES_FIB)
  # Muto-2021-batch-5-small (MES_FIB)
  # Muto-2021-sampled-small (MES_FIB)
  # Muto-2021-small (MES_FIB)
  # Yao-2021-small (NP)
  # paired_root=("HSPC" "OliI" "Basal" "Basal" "Basal" "Basal" "Basal" "Basal")

  pr=HSPC
  for data_name in 10x-Multiome-Pbmc10k-small
  # pr=OliI
  # for data_name in Chen-2019-small
  # pr=Basal
  # for data_name in Ma-2020-batch-53-small Ma-2020-batch-54-small Ma-2020-batch-55-small Ma-2020-batch-56-small Ma-2020-small Ma-2020-sampled-small
  # pr=MES_FIB
  # for data_name in Muto-2021-batch-1-small Muto-2021-batch-2-small Muto-2021-batch-3-small Muto-2021-batch-4-small Muto-2021-batch-5-small
  # pr=NP
  # for data_name in Yao-2021-small

  # for data_name in Ma-2020-batch-55-small Chen-2019-small
  # for data_name in Ma-2020-batch-53-small Ma-2020-batch-54-small Ma-2020-batch-55-small Ma-2020-batch-56-small Ma-2020-small Ma-2020-sampled-small
  # for data_name in 10x-Multiome-Pbmc10k-small Chen-2019-small Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small 
  # for data_name in Ma-2020-small Ma-2020-sampled-small Ma-2020-batch-53-small  Ma-2020-batch-54-small  Ma-2020-batch-55-small  Ma-2020-batch-56-small 
  do
    echo "Running trajectory with $method $data_name"
    CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/trajectory_conservation_metrics_comb.py \
          -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ACTIVE.h5ad \
          -l1 /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-rna.csv \
          -l2 /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name$m-atac.csv\
          -r $pr  \
          -o /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/combo_trajectory_info_$m.yaml
    echo "=================="
  done
done