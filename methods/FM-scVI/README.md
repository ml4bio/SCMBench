## Document of FM-adaptation

To improve the capability of foundation models (FMs) in the specific task of multi-omics integration, we propose a strategy called FM-adaptation. It adds an adaptation module to the end of the foundation model's precomputed embeddings. And in the learning process, only parameters of this adaptation module are updated. In practice, we leverage a basic model, `scVI`, as the adaptation head to test the efficiency of this strategy, and we call this edited model `FM-scVI`. 

We provide scripts in this folder. It requires original measured datasets, precomputed embeddings of each modalities as input.

Example Usage:

```
python run_gpt_scvi.py \
--input-rna RNA.h5ad \
--input-atac ATAC.h5ad \
--rna-pre human-rna.csv \
--atac-pre human-atac.csv \
--output-path output_folder 

# rna-pre and atac-pre: precomputed embedding files
# output-path: output folder path, final embeddings will be generated in this folder.
```
