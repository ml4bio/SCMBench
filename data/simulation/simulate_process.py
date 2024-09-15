import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--path")
args = parser.parse_args()
rna = sc.read_h5ad(args.path+"-RNA.h5ad")
atac = sc.read_h5ad(args.path+"-ATAC.h5ad")
atac.obs.index.name, atac.var.index.name = "cells", "peaks"

sc.pp.filter_genes(rna, min_counts=1)

sc.pp.filter_genes(atac, min_counts=1)

sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

# sc.pp.log1p(rna)
# rna.X = np.array(rna.X,dtype=int)
# temp = rna.X
# temp[temp > 100000] = 100000
# rna.X = temp

# temp = atac.X
# temp[temp > 100000] = 100000
# atac.X = temp
# sc.pp.log1p(atac)
# atac.X = np.array(atac.X,dtype=int)

rna.write_h5ad(args.path +"-RNA.h5ad", compression="gzip")
atac.write_h5ad(args.path +"-ATAC.h5ad", compression="gzip")

# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_2_num_gene_1000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_3_num_gene_1000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_5_num_gene_1000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_2_num_gene_2000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_3_num_gene_2000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_5_num_gene_2000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_2_num_gene_3000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_3_num_gene_3000-noised-effect-5
# python process.py --path ./data/grn_1139_simulated_num_cell_1000_num_batch_5_num_gene_3000-noised-effect-5