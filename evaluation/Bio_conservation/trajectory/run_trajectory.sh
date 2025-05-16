echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++        Integration_trajectory real: '$method'        +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


pr=HSPC

# Root cell type of different datasets:
# PBMC-10x: pr=HSPC
# Chen-2019: pr=OliI
# Ma-2020: pr=Basal
# Muto-2021: pr=MES_FIB
# Yao-2021: pr=NP

for method in bindSC
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




