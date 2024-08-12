# loop over methods on paired datasets
#methods=("bindsc" "Cobolt" "Deepmaps" "GLUE" "LIGER" "MOFA" "Pamona" "PCA" "scJoint" "scMoMaT" "seurat4" "seurat5" "TotalVI" "scVI" "scMDC" "iNMF")
methods=('scJoint' "Cobolt")
HP_setting="default"
#paired_datasets=("10x-Multiome-Pbmc10k-small" "Chen-2019-small" "Ma-2020-batch-53-small" "Ma-2020-batch-54-small" "Ma-2020-batch-55-small" "Ma-2020-batch-56-small" "Ma-2020-sampled-small") #"10x-Multiome-Pbmc10k" "Chen-2019" 
paired_root=("HSPC" "OliI" "Basal" "Basal" "Basal" "Basal" "Basal" "Basal")
unpaired_datasets=("Yao-2021-small")
embedding_path_base="/data2/wangxuesong/benchmark_data/benchmark_embedding/"
unpaired_root=("NP")
data_path_base="/data2/wangxuesong/benchmark_data/new_data/paired_unpaired/new_data/unpaired/"

i=0

for data in "${unpaired_datasets[@]}"; do
    for method in "${methods[@]}"; do
        #for dataset in "${datasets[@]}"; do
            embedding_path="${embedding_path_base}${data}/"
            data_path="${data_path_base}${data}/"
            echo "${method}"
            echo "${data}"
            echo "${unpaired_root[${i}]}"
            python workflow/scripts/trajectory_conservation_metrics.py -d "${data_path}${data}-RNA.h5ad" -l "${embedding_path}${method}_rna.csv" -r "${unpaired_root[${i}]}"  -o "${data}-${method}.yaml"
                    
        #done
    done
    let "i++"
done 