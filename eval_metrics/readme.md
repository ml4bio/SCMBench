### Evaluation metrics 

There are three groups of metrics: biomarker, cell_integration, and trajectory. 

1. Biomarker: 
DARs (differential accessible regions), enriched motifs (identified from DARs), and rna (differential expression genes). 
Two groups of comparison are given for each category: compare_across_methods (comparison between all pairs of methods) and compare with ground truth (calculated from annotated cell type labels). The 'compare_with_prior_markers' in 'rna' group means comparing infered DEGs with existing prior knowledge from biomarker database.

2. Cell integration:
Measures level of cell type matching. We group the evaluation results by different datasets. 

3. Trajectory:
Measures the conservation of infered trajectory from integrated embeddings. We group the evaluation results by different datasets. 