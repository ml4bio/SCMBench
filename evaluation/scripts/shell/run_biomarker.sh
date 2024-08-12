echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Biomarker: scVI         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

method=scVI #scVI iNMF scMDC
# for data_name in Chen-2019 Ma-2020 Yao-2021
for data_name in Ma-2020
# for data_name in 10x-Multiome-Pbmc10k Chen-2019 Ma-2020 Yao-2021
do
  echo "Running exp with $data_name"
  CUDA_VISIBLE_DEVICES=1 python /mnt/nas/user/yixuan/Multiomics-benchmark-main/experiments/BiomarkerDetect/run_biomDetect.py \
        -d /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-RNA.h5ad /mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/$data_name/$data_name-ATAC.h5ad\
        -l /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-rna.csv  /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/$data_name-atac.csv \
        --output-dotplot /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/BiomarkerDetect/figures/output_dotplot.pdf \
        --output-violin /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/BiomarkerDetect/figures/output_violin.pdf  \
        --output-info /mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/$method-output/$data_name/BiomarkerDetect/output_info.yaml
  echo "=================="
done


