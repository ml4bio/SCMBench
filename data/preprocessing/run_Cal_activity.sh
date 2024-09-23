echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '++++++++++++         Multiomics: Activity         +++++++++++++'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '++++++++++++          Activity batch human         ++++++++++++'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

for data_name in 10x-Multiome-Pbmc10k Muto-2021
do
  echo "Running exp with Human $data_name"
  Rscript ./Cal_activity.R \
        --input-rna ../download/$data_name/$data_name-small-RNA.h5ad \
        --input-atac ../download/$data_name/$data_name-small-ATAC.h5ad \
        --output-active ../download/$data_name/$data_name-small-ACTIVE.h5ad  \
        --genome GRCh38 
  echo "=================="
done

# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '++++++++++++        Activity batch Muto-2021       ++++++++++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# for data_name in Muto-2021-sampled Muto-2021-batch-1 Muto-2021-batch-2 Muto-2021-batch-3 Muto-2021-batch-4 Muto-2021-batch-5
# do
#   echo "Running exp with Human $data_name"
#   Rscript ./Cal_activity.R \
#         --input-rna ../download/Muto-2021/$data_name-small-RNA.h5ad \
#         --input-atac ../download/Muto-2021/$data_name-small-ATAC.h5ad \
#         --output-active ../download/Muto-2021/$data_name-small-ACTIVE.h5ad  \
#         --genome GRCh38 
#   echo "=================="
# done

# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '++++++++++++          Activity batch mouse         ++++++++++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# for data_name in Chen-2019 Ma-2020 Yao-2021
# do
#   echo "Running exp with Mouse $data_name"
#   Rscript ./Cal_activity.R \
#         --input-rna ../download/$data_name/$data_name-small-RNA.h5ad \
#         --input-atac ../download/$data_name/$data_name-small-ATAC.h5ad \
#         --output-active ../download/$data_name/$data_name-small-ACTIVE.h5ad  \
#         --genome GRCm38 
#   echo "=================="
# done

# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# echo '++++++++++++         Activity batch Ma-2021        ++++++++++++'
# echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# for data_name in Ma-2020-sampled Ma-2020-batch-53  Ma-2020-batch-54  Ma-2020-batch-55  Ma-2020-batch-56
# do
#   echo "Running exp with Ma-2020 Mouse $data_name"
#   Rscript ./Cal_activity.R \
#         --input-rna ../download/Ma-2020/$data_name-small-RNA.h5ad \
#         --input-atac ../download/Ma-2020/$data_name-small-ATAC.h5ad \
#         --output-active ../download/Ma-2020/$data_name-small-ACTIVE.h5ad  \
#         --genome GRCm38 
#   echo "=================="
# done