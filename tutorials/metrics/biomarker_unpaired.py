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

def jaccard(first, second):
    # Check if 'first' and 'second' are array-like
    if not isinstance(first, (np.ndarray, pd.Series)) or not isinstance(second, (np.ndarray, pd.Series)):
        return 1.0 if first == second else 0.0

    # Convert to set if not entirely null
    first_set = set(first) if not pd.isnull(first).all() else set()
    second_set = set(second) if not pd.isnull(second).all() else set()

    # Compute Jaccard index
    if not first_set and not second_set:
        return 1.0
    return len(first_set.intersection(second_set)) / len(first_set.union(second_set))


parent_dir = '/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download'
latent_dir = '/mnt/nas/user/yixuan/Multiomics-benchmark-main/results'
# datasets = ['Muto-2021-batch-4-small', 'Yao-2021-small', 'Muto-2021-batch-2-small', 'Muto-2021-small']
# datasets = ['Muto-2021-batch-1-small', 'Muto-2021-batch-2-small', 'Muto-2021-batch-3-small', 'Muto-2021-batch-4-small', 'Muto-2021-batch-5-small''Yao-2021-small']
datasets = ['Muto-2021-batch-1-small']
method = [
            'GLUE', 'LIGER','scMoMaT', 'Cobolt','Pamona', 'bindsc',
            'PCA', 'scJoint', 'seurat4','iNMF', 'scVI',
            'scGPT_kidney_zero',
            'UCE','Geneformer','scFoundation'
            ]


for dataset in datasets:
    input_datasets = [f'{parent_dir}/{dataset}/{dataset}-RNA.h5ad', f'{parent_dir}/{dataset}/{dataset}-ATAC.h5ad']


    datasets = [anndata.read_h5ad(item) for item in input_datasets]
    cell_types = [dataset.obs["cell_type"].to_numpy() for dataset in datasets]
    domains = [dataset.obs["domain"].to_numpy() for dataset in datasets]
    vars = [dataset.var_names for dataset in datasets]
    unis = [np.array(dataset.X.todense()) for dataset in datasets]
    cell_types_list = set(list(cell_types[0]))
    print(cell_types_list)
    marker_dict = {}
    
    for main_type in cell_types_list:
        try: 
            for item in method:
                try:
                    input_latents = [f'{latent_dir}/{dataset}/{item}_rna.csv',f'{latent_dir}/{dataset}/{item}_atac.csv']
                    latents = [pd.read_csv(item, header=None, index_col=0).to_numpy() for item in input_latents]
                    masks = [np.apply_along_axis(lambda x: ~np.any(np.isnan(x)), 1, latent) for latent in latents]
                    combined_cell_type = np.concatenate([cell_type[mask] for cell_type, mask in zip(cell_types, masks)])
                    combined_domain = np.concatenate([domain[mask] for domain, mask in zip(domains, masks)])
                    combined_var = np.concatenate([var for var in vars])
                    combined_uni = np.concatenate([uni[mask] for uni, mask in zip(unis, masks)], axis=1)
                    for i, mask in enumerate(masks):
                        rm_pct = 100 * (1 - mask.sum() / mask.size)
                        if rm_pct:
                            print(f"Ignoring {rm_pct:.1f}% cells in dataset {i} due to missing values!")
                    combined_latent = np.concatenate([latent[mask] for latent, mask in zip(latents, masks)], axis=1)
                    combined_predicted_cell_type = classifier_knn(combined_latent, cell_types[0])
                    
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