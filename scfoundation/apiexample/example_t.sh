token=biomape8517a2108634588b383e385e70364d4  ## change to yours
# ### Cell embedding
# taskname=Baron_demo
# tgthighres=a5
# mkdir -p ./demo/${taskname}/${tgthighres}

# python ./client.py --input_type singlecell --output_type cell --pool_type all --pre_normalized F --version 0.2 --tgthighres $tgthighres --data_path ./data/baron_human_samp_19264_fromsaver_demo.csv --save_path ./demo/${taskname}/${tgthighres}/ --token ${token}

# Cell embedding
# taskname=Multiomics_demo
# tgthighres=a5
# # for slice in 1
# for slice in $(seq 1 26)
# do
#     echo "Running exp with $slice"
#     mkdir -p ./demo/${taskname}/act/muto1/${slice}
#     python ./client.py --input_type singlecell --output_type cell --pool_type all --pre_normalized T --version 0.2 --tgthighres $tgthighres --data_path ./data/act/muto1/10x_file_${slice}.csv --save_path ./demo/${taskname}/act/muto1/${slice}/ --token ${token}
# done

## Cell embedding
taskname=Multiomics_demo
tgthighres=a5
for slice in $(seq 1 101)
# for slice in 1 2 3 4 5 6 7 8 9 10
do
    echo "Running exp with $slice"
    mkdir -p ./demo/${taskname}/rna/muto/${slice}
    python ./client.py --input_type singlecell --output_type cell --pool_type all --pre_normalized F --version 0.2 --tgthighres $tgthighres --data_path ./data/rna/muto/10x_file_${slice}.csv --save_path ./demo/${taskname}/rna/muto/${slice}/ --token ${token}
done

## Cell embedding
# taskname=Multiomics_demo
# tgthighres=a5
# mkdir -p ./demo/${taskname}/${tgthighres}
# python ./client.py --input_type singlecell --output_type cell --pool_type all --pre_normalized T --version 0.2 --tgthighres $tgthighres --data_path ./data/10x_act_100.csv --save_path ./demo/${taskname}/${tgthighres}/ --token ${token}





# ### Bulk embedding
# taskname=SCAD_bulk_Etoposide_demo
# tgthighres=f1
# mkdir -p ./demo/${taskname}/${tgthighres}

# python ./client.py --input_type bulk --output_type cell --pool_type all --pre_normalized F --version 0.1 --tgthighres $tgthighres --data_path ./data/Source_exprs_resp_19264.Etoposide_demo.csv --save_path ./demo/${taskname}/${tgthighres}/ --token ${token}

# ### Gene embedding
# taskname=GEARS_demo_batch
# tgthighres=f1
# mkdir -p ./demo/${taskname}/${tgthighres}

# python ./client.py --input_type singlecell --output_type gene --pool_type all --pre_normalized A --version 0.1 --tgthighres $tgthighres --data_path ./data/gene_batch.npy --save_path ./demo/${taskname}/${tgthighres}/ --token ${token}