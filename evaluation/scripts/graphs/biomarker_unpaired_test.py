import os
import random
import warnings
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier

seed_value = 0
os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=anndata._core.views.ImplicitModificationWarning)

def classifier_knn(x, y, neighbor_frac=0.01, test_size=0.8):
    k = max(round(len(y) * neighbor_frac), 1)
    DATA_tr, DATA_te, tr_LABEL, te_LABEL = train_test_split(x, y, test_size=test_size, random_state=42)
    clf = KNeighborsClassifier(n_neighbors=k)
    clf.fit(DATA_tr, tr_LABEL)
    return clf.predict(x)

# def jaccard(first, second):
#     # Check if 'first' and 'second' are array-like
#     if not isinstance(first, (np.ndarray, pd.Series)) or not isinstance(second, (np.ndarray, pd.Series)):
#         return 1.0 if first == second else 0.0

#     # Convert to set if not entirely null
#     first_set = set(first) if not pd.isnull(first).all() else set()
#     second_set = set(second) if not pd.isnull(second).all() else set()

#     # Compute Jaccard index
#     if not first_set and not second_set:
#         return 1.0
#     return len(first_set.intersection(second_set)) / len(first_set.union(second_set))

def jaccard(first, second):
    # Check and convert first input
    if isinstance(first, set):
        first_set = first
    elif isinstance(first, pd.Index):
        first_set = set(first) if not first.empty else set()
    elif isinstance(first, list):
        first_set = set(first) if len(first) else set()
    else:
        first_set = set()

    # Check and convert second input
    if isinstance(second, set):
        second_set = second
    elif isinstance(second, pd.Index):
        second_set = set(second) if not second.empty else set()
    elif isinstance(second, list):
        second_set = set(second) if len(second) else set()
    else:
        second_set = set()
    # Check if both indices are empty
    if not first_set and not second_set:
        # Handle the case where both sets are empty
        # Return 1 or 0 based on your convention for similarity of empty sets
        return 1.0
    return len(first_set.intersection(second_set)) / len(first_set.union(second_set))

parent_dir = '/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download'
latent_dir = '/mnt/nas/user/yixuan/Multiomics-benchmark-main/results'
# datasets = ['Muto-2021-batch-4-small', 'Yao-2021-small', 'Muto-2021-batch-2-small', 'Muto-2021-small']
# datasets = ['Muto-2021-batch-1-small', 'Muto-2021-batch-2-small', 'Muto-2021-batch-3-small', 'Muto-2021-batch-4-small', 'Muto-2021-batch-5-small''Yao-2021-small']
# datasets = ['Muto-2021-batch-1-small']
datasets = ['Yao-2021-small']
method = [
            'GLUE', 'LIGER','scMoMaT', 'Cobolt','Pamona', 'bindsc',
            'PCA', 'scJoint', 'seurat4','iNMF', 'scVI',
            'scGPT',
            'UCE','Geneformer','scFoundation'
            ]


for dataset in datasets:
    input_datasets = f'{parent_dir}/{dataset}/{dataset}-RNA.h5ad'
    data = anndata.read_h5ad(input_datasets)
    cell_type, domain, var = data.obs["cell_type"].to_numpy(), data.obs["domain"].to_numpy(), data.var_names
    uni = np.array(data.X.todense())
    cell_types_list = set(cell_type)
    print(cell_types_list)
    marker_dict = {}
    # cell_types_list=['ENDO']
    for main_type in cell_types_list:
        try: 
            for item in method:
                try:
                    input_latent = f'{latent_dir}/{dataset}/{item}_rna.csv'
                    latent = pd.read_csv(input_latent, header=None, index_col=0).to_numpy()
                    mask = ~np.isnan(latent).any(axis=1)
                    combined_latent = latent[mask]
                    predicted_cell_type = classifier_knn(combined_latent, cell_type[mask])

                    adata = sc.AnnData(X=uni[mask], obs={'domain': domain[mask], 'cell_type': cell_type[mask], 'predicted_cell_type': predicted_cell_type}, var=var.to_frame(index=False))
                    sc.pp.normalize_total(adata)
                    sc.pp.log1p(adata)
                    sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="seurat")
                    adata = adata[:, adata.var.highly_variable]

                    adata.obs['predicted_cell_type'] = np.array(predicted_cell_type)

                    if main_type in adata.obs['predicted_cell_type'].unique():
                        sc.tl.rank_genes_groups(adata, groupby='predicted_cell_type', groups=[main_type], method='t-test', n_genes=100)
                        marker_dict[item] = sc.get.rank_genes_groups_df(adata, group=main_type).set_index('names').index.tolist()

                    else:
                        print(f"Error: '{main_type}' not present in 'predicted_cell_type'")
                        marker_dict[item] = set()
                    print("Finish:", item)
                except (FileNotFoundError, ValueError, IndexError) as e:
                    print(f"{e.__class__.__name__} for {item}: {str(e)}")
                    continue
            plot_directory = f'/mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/graphs/{dataset}'
            os.makedirs(plot_directory, exist_ok=True)
            np.save(f'{plot_directory}/marker_dict_{main_type}.npy', marker_dict)

            result_dict = {(k, l): jaccard(marker_dict[k], marker_dict[l]) for k in marker_dict.keys() for l in marker_dict.keys()}
            data = np.array(list(result_dict.values())).reshape(len(marker_dict.keys()), len(marker_dict.keys()))
            np.save(f'{plot_directory}/np_jaccard_{main_type}.npy', data)

            medians = np.median(data, axis=1)
            np.save(f'{plot_directory}/np_median_{main_type}.npy', medians)

            sorted_indices = np.argsort(medians)[::-1]
            sorted_medians = medians[sorted_indices]
            sorted_methods = np.array(list(marker_dict.keys()))[sorted_indices]

            fig, ax = plt.subplots(figsize=(9, 9))
            mask = np.tril(np.ones_like(data, dtype=bool))
            sns.heatmap(data, mask=mask, annot=True, cmap=sns.diverging_palette(220, 20, as_cmap=True), ax=ax, cbar=False)
            ax.xaxis.tick_top()
            ax.set_xticklabels(marker_dict.keys(), rotation=45, fontsize=14)
            ax.yaxis.tick_right()
            ax.set_yticklabels(marker_dict.keys(), rotation=0, fontsize=14)
            plt.tight_layout()
            plt.savefig(f'{plot_directory}/{main_type}_heatmap.pdf', transparent=True)

            fig, ax = plt.subplots(figsize=(7, 4))
            ax.bar(np.arange(len(sorted_methods)), sorted_medians, width=0.4, color=[sns.diverging_palette(220, 20)[4] if m >= 0.9 else (sns.diverging_palette(220, 20)[1] if m >= 0.8 else sns.diverging_palette(220, 20)[0]) for m in sorted_medians])
            ax.set_xticks(np.arange(len(sorted_methods)))
            ax.set_xticklabels(sorted_methods, rotation=60, fontsize=12)
            ax.set_ylim(0, 1)
            ax.set_ylabel('Median', fontsize=16)
            plt.tight_layout()
        except (FileNotFoundError, ValueError, IndexError) as e:
            print(f"{e.__class__.__name__} for {item}: {str(e)}")
            continue