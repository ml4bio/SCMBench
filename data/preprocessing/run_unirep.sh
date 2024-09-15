echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '++++++++++++++++++         unirep         ++++++++++++++++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

for data_name in 10x-Multiome-Pbmc10k Chen-2019 Ma-2020 Muto-2021 Yao-2021
do
      echo "Running exp with $data_name -ATAC"
      CUDA_VISIBLE_DEVICES=1 python ./atac_unirep.py \
            -i ../download/$data_name/$data_name-small-ATAC.h5ad \
            -o ../download/$data_name/$data_name-small-ATAC_uni.h5ad 
      echo "=================="

      echo "Running exp with $data_name -RNA"
      CUDA_VISIBLE_DEVICES=1 python ./rna_unirep.py \
            -i ../download/$data_name/$data_name-small-RNA.h5ad\
            -o ../download/$data_name/$data_name-small-RNA_uni.h5ad
      echo "==================="
done

echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '++++++++++++         unirep batch Ma-2020        ++++++++++++'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

for data_name in Ma-2020-sampled Ma-2020-batch-53  Ma-2020-batch-54  Ma-2020-batch-55  Ma-2020-batch-56 
do
      echo "Running exp with $data_name -ATAC"
      CUDA_VISIBLE_DEVICES=1 python ./atac_unirep.py \
            -i ../download/Ma-2020/$data_name-small-ATAC.h5ad \
            -o ../download/Ma-2020/$data_name-small-ATAC_uni.h5ad 
      echo "=================="

      echo "Running exp with $data_name -RNA"
      CUDA_VISIBLE_DEVICES=1 python ./rna_unirep.py \
            -i ../download/Ma-2020/$data_name-small-RNA.h5ad\
            -o ../download/Ma-2020/$data_name-small-RNA_uni.h5ad
      echo "==================="
done


echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo '++++++++++++         unirep batch Muto-2021        ++++++++++++'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

for data_name in Muto-2021-sampled Muto-2021-batch-1 Muto-2021-batch-2 Muto-2021-batch-3 Muto-2021-batch-4 Muto-2021-batch-5 
do
      echo "Running exp with $data_name -ATAC"
      CUDA_VISIBLE_DEVICES=1 python ./atac_unirep.py \
            -i ../download/Muto-2021/$data_name-small-ATAC.h5ad \
            -o ../download/Muto-2021/$data_name-small-ATAC_uni.h5ad 
      echo "=================="

      echo "Running exp with $data_name -RNA"
      CUDA_VISIBLE_DEVICES=1 python ./rna_unirep.py \
            -i ../download/Muto-2021/$data_name-small-RNA.h5ad\
            -o ../download/Muto-2021/$data_name-small-RNA_uni.h5ad
      echo "==================="
done
