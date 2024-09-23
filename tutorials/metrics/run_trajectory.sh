echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++        Integration_trajectory real: '$method'        +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
pr=HSPC
# for method in bindSC Cobolt Deepmaps Geneformer GLUE iNMF LIGER MMD_MA MOFA Pamona PCA scFoundation 
for method in scGPT-bc scGPT-bc-zero scGPT scGPT-zero scJoint scMDC scMoMaT scVI Seurat4 Seurat5 TotalVI UCE UnionCom
# for method in bindSC
do
  for data_name in 10x-Multiome-Pbmc10k
  do
    echo "Running trajectory with $method $data_name"
    python trajectory_conservation_joint.py \
          -d ../../data/download/$data_name/$data_name-small-RNA.h5ad ../../data/download/$data_name/$data_name-small-ACTIVE.h5ad \
          -l ../../results/$method-output/$data_name/$data_name-small-rna.csv  ../../results/$method-output/$data_name/$data_name-small-atac.csv \
          -r $pr  \
          -c ../analysis/traj/raw_combine_traj.h5ad \
          -o ../analysis/traj/$method-combo_trajectory_info_active_all.yaml
    echo "=================="
  done
done


# pr=HSPC
# # for method in bindSC Cobolt Deepmaps Geneformer GLUE iNMF LIGER MMD_MA MOFA Pamona PCA scFoundation scGPT-bc scGPT-bc-zero scGPT scGPT-zero scJoint scMDC scMoMaT scVI Seurat4 Seurat5 TotalVI UCE UnionCom
# for method in bindSC
# do
#   for data_name in 10x-Multiome-Pbmc10k
#   do
#     echo "Running trajectory with $method $data_name"
#     python trajectory_conservation_joint.py \
#           -d ../../data/download/$data_name/$data_name-small-RNA.h5ad ../../data/download/$data_name/$data_name-small-ACTIVE.h5ad \
#           -l ../../results/$method-output/$data_name/$data_name-small-rna.csv  ../../results/$method-output/$data_name/$data_name-small-atac.csv \
#           -r $pr  \
#           -b \
#           -bn 'HSPC' 'CD8 TEM_2' 'CD8 TEM_1' 'CD4 Naive' 'CD8 Naive' 'CD4 TCM' 'CD4 TEM' 'Naive B' 'Memory B' 'Intermediate B' \
#           -c ../analysis/traj/raw_combine_traj_branch1_idx4.h5ad \
#           -o ../analysis/traj/$method-combo_trajectory_info_active_branch.yaml
#     echo "=================="
#   done
# done

