## Evaluation Section

Codes for downstream evaluation tasks, including Multi-Omic Integration Accuracy, Bio-Conservation, and Batch effect & detecting over-correction. We provide the python script for each metrics in this folder for evaluation.

### Multi-Omic Integration Accuracy
Run [Integration_accuracy/integration_accuracy.py](Integration_accuracy/integration_accuracy.py) with shell script [Integration_accuracy/run_integration.sh](Integration_accuracy/run_integration.sh) for multi-omic integration evaluation. The input of integration accuracy script are `RNA_uni.adata`, `ATAC_uni.adata`, `RNA_emb.csv`, `ATAC_emb.csv`. `RNA_uni.adata` and `ATAC_uni.adata` are unimodal representations calculated from [../data/preprocessing/run_unirep.sh](../data/preprocessing/run_unirep.sh). 

Example Usage:

```
python integration_accuracy.py \
-d RNA_uni.h5ad ATAC_uni.h5ad \ 
-l rna.csv atac.csv \ 
-o cell_integration_info.yaml

# -d: Unimodel representations of RNA and ATAC dataset
# -l: Inferred latent embeddings for rna and atac
```

Output contins:
- Mean Average Precision (MAP)
- Normalized Mutual Info (NMI)
- Average Silhouette Width (ASW)
- Adjusted Rand Index (ARI)

### Bio-Conservation

#### Biomarker and DARs

Run [Bio_conservation/biomarker/biomarker.py](Bio_conservation/biomarker/biomarker.py) for measuring the conservation score of biological features, marker genes, and differential accessible regions (DARs). The biomarker conservation score is evaluated across multiple methods, so you need to input multiple embeddings at once. The script will calculate the Jaccard Similarity Index (JSI) for each pair of methods. \

Example Usage:

```
python biomarker/biomarker.py
        -d RNA.h5ad \ # Input RNA dataset, for DARs evaluation use ATAC dataset
        -l Embeddings \ # Path of Folder containing embeddings.csv
        -m scVI GLUE LIGER \ # List of methods you want to evaluate
        --mode rna \ # 'rna' for marker genes, 'atac' for DARs
        -o output_folder # Path of output folder
```

Input file format:
```
Embeddings
├── scVI_rna.csv    
├── scVI_atac.csv     
├── GLUE_rna.csv   
├── GLUE_atac.csv   
├── LIGER_rna.csv    
├── LIGER_atac.csv    
...
```

Output:
- Marker gene or DARs dictionary and JSI matrix (shape n*n for n methods) for each cell type
- Conservation scores of all input methods

#### Enriched Motifs

Run [Bio_conservation/motifs/run_gimmemotifs.sh](Bio_conservation/motifs/run_gimmemotifs.sh)

#### Trajectory conservation
Run [Bio_conservation/trajectory/trajectory_conservation_joint.py](Bio_conservation/trajectory/trajectory_conservation_joint.py) with [Bio_conservation/trajectory/run_trajectory.sh](Bio_conservation/trajectory/run_trajectory.sh) for trajectory conservation. It requires `ACTIVE.adata` as input, which is the gene activity matrix calculated from `ATAC.adata` from [../data/preprocessing/Cal_activity.R](../data/preprocessing/Cal_activity.R). 

```
python trajectory/trajectory_conservation_joint.py \
        -d RNA.h5ad ACTIVE.h5ad \ 
        -l rna.csv  atac.csv \ # inferred latent embeddings
        -r HSPC \ # root cell type, set as HSPC for PBMC-10x
        -c raw_combine_traj.h5ad \ # output conbined anndata
        -o combo_trajectory_info_active_all.yaml 
```

Output:
- Trajectory Conservation Score.

### Batch Correction
Run [Batch_correction/batch_effect_correction.py](Batch_correction/batch_effect_correction.py) for batch correction evaluation. 

Example Usage:

```
python Batch_correction/batch_effect_metrics.py \
    -d RNA_uni.h5ad ATAC_uni.h5ad \
    -l rna.csv atac.csv \
    -o batch_effect_info.yaml 
```

Output:
- Batch-Average Silhouette Width (Batch-ASW)
- iLISI
- kBET
- Graph Connectivity
