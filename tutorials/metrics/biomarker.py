import scanpy as sc
import seaborn as sns
import os 
from matplotlib import rcParams
import matplotlib.pyplot as plt
import pickle
import argparse
import functools
import pathlib
import warnings
from itertools import compress

import anndata
import numpy as np
import pandas as pd
import yaml

from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.spatial
import sklearn.metrics
import sklearn.neighbors
from anndata import AnnData
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

def classifer_knn(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01, test_size: float = 0.8, **kwargs
) -> float:
    k = max(round(y.shape[0] * neighbor_frac), 1)
    DATA_tr, DATA_te, tr_LABEL, te_LABEL = train_test_split(x, y, test_size=test_size, random_state=42)
    clf = KNeighborsClassifier(n_neighbors=k)
    clf.fit(DATA_tr, tr_LABEL)
    L_pred = clf.predict(x)
    return L_pred

# def jaccard(first, second):
#     return len(set(first).intersection(second)) / len(set(first).union(second))
def jaccard(first, second):
    # Check and convert first input
    if isinstance(first, set):
        first_set = first
    elif isinstance(first, pd.Index):
        first_set = set(first) if not first.empty else set()
    elif isinstance(first,np.ndarray):
        first_set = set(first)
    else:
        first_set = set()

    # Check and convert second input
    if isinstance(second, set):
        second_set = second
    elif isinstance(second, pd.Index):
        second_set = set(second) if not second.empty else set()
    elif isinstance(second,np.ndarray):
        second_set = set(second)
    else:
        second_set = set()
    # Check if both indices are empty
    if not first_set and not second_set:
        # Handle the case where both sets are empty
        # Return 1 or 0 based on your convention for similarity of empty sets
        return 1.0
    return len(first_set.intersection(second_set)) / len(first_set.union(second_set))

def precision(first, second):
    # print(type(first))    
    # Check and convert first input
    if isinstance(first, set):
        first_set = first
    elif isinstance(first, pd.Index):
        first_set = set(first) if not first.empty else set()
    else:
        first_set = set()

    # Check and convert second input
    if isinstance(second, set):
        second_set = second
    elif isinstance(second, pd.Index):
        second_set = set(second) if not second.empty else set()
    elif isinstance(second,np.ndarray):
        second_set = set(second)
    else:
        second_set = set()
    # Check if both indices are empty
    if not first_set and not second_set:
        # print('**')
        # Handle the case where both sets are empty
        # Return 1 or 0 based on your convention for similarity of empty sets
        return 1.0
    # print('##')
    return len(first_set.intersection(second_set)) / len(second_set)

def parse_args():
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Compute integration metrics for paired samples"
    )
    parser.add_argument(
        "-d", "--datasets", dest="datasets", type=pathlib.Path, required=True,
        nargs="+", help="Path to datasets (.h5ad)"
    )
    parser.add_argument(
        "-l", "--latent_dir", dest="latent_dir", type=pathlib.Path, required=True,
        help="Fold Directory to latent embeddings"
    )
    parser.add_argument(
        "-m", "--methods", dest="methods", type=str, required=True,
        nargs="+", help="List of evaluated methods"
    )
    parser.add_argument(
        "--mode", dest="mode", type=str, default="rna",
        help="Whether evaluate biomarker (rna) or differential accessible regions (atac)"
    )
    parser.add_argument(
        "--cell-type", dest="cell_type", type=str, default="cell_type",
        help="Column name in obs specifying cell types"
    )
    parser.add_argument(
        "--domain", dest="domain", type=str, default="domain",
        help="Column name in obs specifying domain"
    )
    parser.add_argument(
        "-p", "--paired", dest="paired", action="store_true",
        help="Whether the latent embeddings are paired"
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", type=pathlib.Path, required=True,
        help="Path to output directory"
    )
    return parser.parse_args()

def main(args):
    # required_data = f'Ma-2020-batch-53-small'

    # method = [
    #             'GLUE', 'LIGER', 'TotalVI', 'UnionCom',
    #             'scMoMaT', 'MMD_MA', 'Deepmaps',
    #             'Cobolt', 'scMDC','bindsc',
    #             'Pamona', 'PCA', 'scJoint', 'seurat4',
    #             'seurat5',
    #             'iNMF', 'scVI', 'MOFA','Harmony',
    #             'UCE','scGPT','scGPT_scvi'
    #             ]
    method = args.methods
    input_datasets = args.datasets
    datasets = [anndata.read_h5ad(item) for item in input_datasets]
    cell_types = [dataset.obs[args.cell_type].to_numpy() for dataset in datasets]
    domains = [dataset.obs[args.domain].to_numpy() for dataset in datasets]
    vars = [dataset.var_names for dataset in datasets]
    unis = [np.array(dataset.X.todense()) for dataset in datasets]
    cell_types_list = set(list(cell_types[0]))
    print(cell_types_list)

    bio_scores = []
    for main_type in cell_types_list:
        marker_dict = {}
        adata = datasets[0]
        if len(adata[adata.obs[args.cell_type]==main_type])<5:
            print(f"Error: '{main_type}' not present in 'cell_type'")
            continue
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=5000,
            flavor="seurat",
            batch_key=None,
        )
        adata = adata[:, adata.var.highly_variable]

        for item in method:
            # Considering some methods generate separated embeddings for diverse modalities, rather than joint embedding, 
            # we use the single embedding of corresponding modality for evaluation.
            # For methods generating joint embedding, we take the joint one as input.
            input_latents = [f'{args.latent_dir}/{item}_{args.mode}.csv']
            latents = [pd.read_csv(item, header=None, index_col=0).to_numpy() for item in input_latents]
            masks = [np.apply_along_axis(lambda x: ~np.any(np.isnan(x)), 1, latent) for latent in latents]
            for i, mask in enumerate(masks):
                rm_pct = 100 * (1 - mask.sum() / mask.size)
                if rm_pct:
                    print(f"Ignoring {rm_pct:.1f}% cells in dataset {i} due to missing values!")
            combined_cell_type = np.concatenate([cell_type[mask] for cell_type, mask in zip(cell_types, masks)])
            combined_domain = np.concatenate([domain[mask] for domain, mask in zip(domains, masks)])
            combined_var = np.concatenate([var for var in vars])
            combined_uni = np.concatenate([uni[mask] for uni, mask in zip(unis, masks)], axis=1)
            combined_latent = np.concatenate([latent[mask] for latent, mask in zip(latents, masks)], axis=1)

            # Predict cell types with infered embeddings using KNN
            combined_predicted_cell_type = classifer_knn(combined_latent, cell_types[0])
            adata = sc.AnnData(X=combined_uni,
                                obs=list(domains[0]),
                                var=list(combined_var))
            adata.var_names = combined_var
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=5000,
                flavor="seurat",
                batch_key=None,
            )
            adata = adata[:, adata.var.highly_variable]
            adata.obs['cell_type'] = np.array(cell_types[0])
            adata.obs['predicted_cell_type'] = np.array(combined_predicted_cell_type)

            if len(adata[adata.obs['predicted_cell_type']==main_type])<5:
                print(f"Error: '{main_type}' not present in 'predicted_cell_type'")
                marker_dict[item] = set()
                continue

            # Calculate differential genes or accessible regions
            sc.tl.rank_genes_groups(
                adata, 
                groupby = 'predicted_cell_type', 
                groups = [main_type],
                method = 't-test',
                n_genes = 100,
                )
            md = ( 
                sc.get.rank_genes_groups_df(adata,group=[main_type]) 
                .set_index('names', drop=True) 
            ) 
            marker_dict[item] = md.index
            print("Finish: ", item)

        safe_main_type = main_type.replace('/', '_')

        keys = list(marker_dict.keys())
        result_dict = {}

        # Calculate JSI score between each pair of methods
        for k in keys:
            for l in keys:
                result_dict[(k,l)] = marker_dict.get((l,k), jaccard(marker_dict[k], marker_dict[l]))

        data = np.array(list(result_dict.values()) ).reshape(len(keys),len(keys))

        os.makedirs(args.output_dir, exist_ok=True)

        # np.save(f'{args.output_dir}/{safe_main_type}_jsi.npy',data)
        # with open(f'{args.output_dir}/{safe_main_type}_marker_dict.pkl','wb') as f:
        #     pickle.dump(marker_dict,f)

        # Assuming 'data' and 'keys' variables are already defined
        mask = np.tril(np.ones_like(data, dtype=bool))
        medians = np.median(data, axis=1)  # Calculate medians across all methods
        means = np.mean(data,axis=1) # Calculate mean across all methods
        bio_scores.append(medians)

    # Calculate biomarker scores for all methods by averaging JSI across all cell types and methods
    bio_scores = np.stack(bio_scores,axis=1)
    bio_scores_mean = bio_scores.mean(1)

    bio_dict = dict(zip(keys,bio_scores_mean))
    for k,v in bio_dict.items():
        bio_dict[k] = v.item()

    if args.mode == 'rna':
        eval_mode = 'biomarker'
    if args.mode == 'atac':
        eval_mode = 'DARs'

    # Save results
    with open(f"{args.output}/{eval_mode}_scores.yaml","w") as f:
        yaml.dump(bio_dict, f)

        
if __name__ == '__main__':
    main(parse_args())