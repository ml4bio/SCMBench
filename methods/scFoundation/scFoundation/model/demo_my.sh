# Copyright 2023 BioMap (Beijing) Intelligence Technology Limited

# Enhancement
python get_embedding.py --task_name Bcell --input_type singlecell --output_type cell --pool_type all --tgthighres a5 --data_path /mnt/nas/user/yixuan/Bio_clock/blood/Bcell_pseudo_100/Bcell_19264_4.csv --save_path /mnt/nas/user/yixuan/scFoundation/apiexample/data --pre_normalized T --version rde
# python get_embedding.py --task_name Bcell --input_type singlecell --output_type cell --pool_type all --tgthighres a5 --data_path /mnt/nas/user/yixuan/Bio_clock/blood/Bcell_pseudo_100/Bcell_19264_1.csv --save_path /mnt/nas/user/yixuan/Bio_clock/results/embeddings --pre_normalized F --version rde
# python get_embedding.py --task_name Bcell --input_type bulk --output_type cell --pool_type all --tgthighres f1 --data_path /mnt/nas/user/yixuan/Bio_clock/blood/Bcell_pseudo_100/Bcell_19264_1.csv --save_path /mnt/nas/user/yixuan/Bio_clock/results/embeddings --pre_normalized F --version ce
# python get_embedding.py --task_name SCAD_bulk_PLX4720_451Lu --input_type bulk --output_type cell --pool_type all --tgthighres f1 --data_path ./examples/SCAD/Source_exprs_resp_19264.PLX4720_451Lu.csv --save_path ./examples/SCAD/ --pre_normalized F --version ce --demo

## different foldchange change tgthighres to f1, f1.5, f2 ...
# ! python get_embedding.py --task_name Baron --input_type singlecell --output_type cell --pool_type all --tgthighres f1 --data_path ./examples/Baron_enhancement.csv --save_path ./examples/ --pre_normalized F --version rde
# ! python get_embedding.py --task_name Baron --input_type singlecell --output_type cell --pool_type all --tgthighres f1.5 --data_path ./examples/Baron_enhancement.csv --save_path ./examples/ --pre_normalized F --version rde
# ! python get_embedding.py --task_name Baron --input_type singlecell --output_type cell --pool_type all --tgthighres f2 --data_path ./examples/Baron_enhancement.csv --save_path ./examples/ --pre_normalized F --version rde

# python get_embedding.py --task_name Zheng68K --input_type singlecell --output_type cell --pool_type all --tgthighres f1 --data_path ./examples/enhancement/pbmc68ksorted_count.npz --save_path ./examples/enhancement/ --pre_normalized F --version rde --demo

