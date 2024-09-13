for method in scMDC
do 
  echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  echo '+++++         Integration_trajectory real: '$method'         +++++'
  echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  pr=HSPC
  for data_name in 10x-Multiome-Pbmc10k-small
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