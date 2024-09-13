## Evaluation Section

Codes for downstream evaluation tasks, including Multi-Omic Integration Accuracy, Bio-Conservation, and Batch effect & detecting over-correction. We provide the python script for each metrics in `scripts/` and the corresponding example shell script for running in `scripts/shell`.

Input of all evaluation scripts are `RNA_uni.adata`, `ATAC_uni.adata`, `RNA_emb.csv`, `ATAC_emb.csv`. `RNA_uni.adata` and `ATAC_uni.adata` are unimodal representations calculated from [scripts/shell/run_reginf.sh](scripts/shell/run_reginf.sh).

### Multi-Omic Integration Accuracy
Run [scripts/cell_integration.py](scripts/shell/run_integration.sh) with shell script [scripts/shell/run_integration.sh](scripts/shell/run_integration.sh) for multi-omic integration evaluation.

Output:
- Mean Average Precision (MAP)
- Normalized Mutual Info (NMI)
- Average Silhouette Width (ASW)
- Adjusted Rand Index (ARI)

### Bio-Conservation
Run [/scripts/graphs/biomarker_paired.py](/scripts/graphs/biomarker_paired.py) for paired datasets.\
Run [/scripts/graphs/biomarker_unpaired.py](/scripts/graphs/biomarker_unpaired.py) for unpaired datasets.\
Run [/scripts/trajectory_conservation_metrics](/scripts/trajectory_convervation_metrics.py) for trajectory conservation.

Output:
- Jaccard Similarity Index (JSI) between outputs of two tools.
- Trajectory Conservation Score.

### Batch Correction
Run [scripts/batch_effect_metrics.py](scripts/batch_effect_metrics.py) for batch correction evaluation.

Output:
- Batch-Average Silhouette Width (Batch-ASW)
- iLISI
- kBET
- Graph Connectivity
