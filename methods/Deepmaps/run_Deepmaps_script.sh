echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Multiomics: Deepmaps         +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'

CUDA_VISIBLE_DEVICES=1 /workspace/wangyixuan/.conda/envs/zgy_r/bin/Rscript /mnt/nas/user/yixuan/zgy/comp_effic/Deepmaps/run_Deepmaps.R \
    --input-rna /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k-small/10x-Multiome-Pbmc10k-small-RNA_5000.h5ad \
    --input-atac /mnt/nas/user/yixuan/SCMBench/data/download/10x-Multiome-Pbmc10k-small/10x-Multiome-Pbmc10k-small-ATAC_50000.h5ad \
    --output-rna /mnt/nas/user/yixuan/zgy/comp_effic/Deepmaps/rna_latent.csv \
    --output-atac /mnt/nas/user/yixuan/zgy/comp_effic/Deepmaps/atac_latent.csv \
    --run-info /mnt/nas/user/yixuan/zgy/comp_effic/Deepmaps/run_info_gpu.json \
    --genome GRCh38
