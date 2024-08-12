import logging

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse

class NeighborsError(Exception):
    def __init__(self, message):
        self.message = message
def diffusion_conn(adata, min_k=50, copy=True, max_iterations=26):
    """
    Diffusion for connectivites matrix extension
    This function performs graph diffusion on the connectivities matrix until a
    minimum number `min_k` of entries per row are non-zero.

    Note:
    Due to self-loops min_k-1 non-zero connectivies entries is actually the stopping
    criterion. This is equivalent to `sc.pp.neighbors`.

    Returns:
       The diffusion-enhanced connectivities matrix of a copy of the AnnData object
       with the diffusion-enhanced connectivities matrix is in
       `adata.uns["neighbors"]["conectivities"]`
    """
    if "neighbors" not in adata.uns:
        raise ValueError(
            "`neighbors` not in adata object. " "Please compute a neighbourhood graph!"
        )

    if "connectivities" not in adata.obsp:
        raise ValueError(
            "`connectivities` not in `adata.obsp`. "
            "Please pass an object with connectivities computed!"
        )

    T = adata.obsp["connectivities"]

    # Normalize T with max row sum
    # Note: This keeps the matrix symmetric and ensures |M| doesn't keep growing
    T = scipy.sparse.diags(1 / np.array([T.sum(1).max()] * T.shape[0])) * T

    M = T

    # Check for disconnected component
    n_comp, labs = scipy.sparse.csgraph.connected_components(
        adata.obsp["connectivities"], connection="strong"
    )

    if n_comp > 1:
        tab = pd.value_counts(labs)
        small_comps = tab.index[tab < min_k]
        large_comp_mask = np.array(~pd.Series(labs).isin(small_comps))
    else:
        large_comp_mask = np.array([True] * M.shape[0])

    T_agg = T
    i = 2
    while ((M[large_comp_mask, :][:, large_comp_mask] > 0).sum(1).min() < min_k) and (
        i < max_iterations
    ):
        print(f"Adding diffusion to step {i}")
        T_agg *= T
        M += T_agg
        i += 1

    if (M[large_comp_mask, :][:, large_comp_mask] > 0).sum(1).min() < min_k:
        raise ValueError(
            "could not create diffusion connectivities matrix"
            f"with at least {min_k} non-zero entries in"
            f"{max_iterations} iterations.\n Please increase the"
            "value of max_iterations or reduce k_min.\n"
        )

    M.setdiag(0)

    if copy:
        adata_tmp = adata.copy()
        adata_tmp.uns["neighbors"].update({"diffusion_connectivities": M})
        return adata_tmp

    else:
        return M


def diffusion_nn(adata, k, max_iterations=26):
    """
    Diffusion neighbourhood score
    This function generates a nearest neighbour list from a connectivities matrix
    as supplied by BBKNN or Conos. This allows us to select a consistent number
    of nearest neighbours across all methods.

    Return:
       `k_indices` a numpy.ndarray of the indices of the k-nearest neighbors.
    """
    if "neighbors" not in adata.uns:
        raise ValueError(
            "`neighbors` not in adata object. " "Please compute a neighbourhood graph!"
        )

    if "connectivities" not in adata.obsp:
        raise ValueError(
            "`connectivities` not in `adata.obsp`. "
            "Please pass an object with connectivities computed!"
        )

    T = adata.obsp["connectivities"]

    # Row-normalize T
    T = scipy.sparse.diags(1 / T.sum(1).A.ravel()) * T

    T_agg = T**3
    M = T + T**2 + T_agg
    i = 4

    while ((M > 0).sum(1).min() < (k + 1)) and (i < max_iterations):
        # note: k+1 is used as diag is non-zero (self-loops)
        print(f"Adding diffusion to step {i}")
        T_agg *= T
        M += T_agg
        i += 1

    if (M > 0).sum(1).min() < (k + 1):
        raise NeighborsError(
            f"could not find {k} nearest neighbors in {max_iterations}"
            "diffusion steps.\n Please increase max_iterations or reduce"
            " k.\n"
        )

    M.setdiag(0)
    k_indices = np.argpartition(M.A, -k, axis=1)[:, -k:]

    return k_indices
def kBET(
    adata,
    batch_key,
    label_key,
    type_,
    embed=None,
    scaled=True,
    return_df=False,
    verbose=False,
):
    """kBET score

    Compute the average of k-nearest neighbour batch effect test (kBET) score per label.
    This is a wrapper function of the implementation by `Büttner et al. 2019`_.
    kBET measures the bias of a batch variable in the kNN graph.
    Specifically, kBET is quantified as the average rejection rate of Chi-squared tests of local vs  global batch label
    distributions.
    This means that smaller values indicate better batch mixing.
    By default the original kBET score is scaled between 0 and 1 so that larger scores are associated with better batch
    mixing.

    .. _Büttner et al. 2019: https://doi.org/10.1038/s41592-018-0254-1

    :param adata: anndata object to compute kBET on
    :param batch_key: name of batch column in adata.obs
    :param label_key: name of cell identity labels column in adata.obs
    :param `type_`: type of data integration, one of 'knn', 'embed' or 'full'
    :param embed: embedding key in ``adata.obsm`` for embedding and feature input
    :param scaled: whether to scale between 0 and 1
        with 0 meaning low batch mixing and 1 meaning optimal batch mixing
        if scaled=False, 0 means optimal batch mixing and 1 means low batch mixing
    :return:
        kBET score (average of kBET per label) based on observed rejection rate.
        If ``return_df=True``, also return a ``pd.DataFrame`` with kBET observed
        rejection rate per cluster

    This function can be applied to all integration output types and recomputes the kNN graph for feature and embedding
    output with specific parameters.
    Thus, no preprocessing is required, but the correct output type must be specified in ``type_``.

    **Examples**

    .. code-block:: python

        # full feature integration output or unintegrated data
        scib.me.kBET(
            adata, batch_key="batch", label_key="celltype", type_="full", embed="X_pca"
        )

        # embedding output
        scib.me.kBET(
            adata, batch_key="batch", label_key="celltype", type_="embed", embed="X_emb"
        )

        # kNN output
        scib.me.kBET(adata, batch_key="batch", label_key="celltype", type_="knn")

    """
    try:
        import rpy2.rinterface_lib.callbacks
        import rpy2.rinterface_lib.embedded
        import rpy2.robjects as ro

        rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError("rpy2 not installed, do: pip install rpy2")



    try:
        ro.r("library(kBET)")
    except rpy2.rinterface_lib.embedded.RRuntimeError as e:
        raise ModuleNotFoundError("kBET not installed, do: devtools::install_github('theislab/kBET')")

    # compute connectivities for non-knn type data integrations
    # and increase neighborhoods for knn type data integrations
    if type_ != "knn" and embed is not None:
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=50, use_rep=embed, copy=True)
    else:
        # check if pre-computed neighbours are stored in input file
        adata_tmp = adata.copy()
        if "diffusion_connectivities" not in adata.uns["neighbors"]:
            if verbose:
                print("Compute diffusion neighbours")
            adata_tmp = diffusion_conn(adata, min_k=50, copy=True)
        adata_tmp.obsp["connectivities"] = adata_tmp.uns["neighbors"][
            "diffusion_connectivities"
        ]

    if verbose:
        print(f"batch: {batch_key}")

    # set upper bound for k0
    size_max = 2**31 - 1

    # prepare call of kBET per cluster
    kBET_scores = {"cluster": [], "kBET": []}
    for clus in adata_tmp.obs[label_key].unique():

        # subset by label
        adata_sub = adata_tmp[adata_tmp.obs[label_key] == clus, :].copy()

        # check if neighborhood size too small or only one batch in subset
        if np.logical_or(
            adata_sub.n_obs < 10, len(adata_sub.obs[batch_key].cat.categories) == 1
        ):
            print(f"{clus} consists of a single batch or is too small. Skip.")
            score = np.nan
        else:
            quarter_mean = np.floor(
                np.mean(adata_sub.obs[batch_key].value_counts()) / 4
            ).astype("int")
            k0 = np.min([70, np.max([10, quarter_mean])])
            # check k0 for reasonability
            if k0 * adata_sub.n_obs >= size_max:
                k0 = np.floor(size_max / adata_sub.n_obs).astype("int")

            matrix = np.zeros(shape=(adata_sub.n_obs, k0 + 1))

            if verbose:
                print(f"Use {k0} nearest neighbors.")
            n_comp, labs = scipy.sparse.csgraph.connected_components(
                adata_sub.obsp["connectivities"], connection="strong"
            )

            if n_comp == 1:  # a single component to compute kBET on
                try:
                    nn_index_tmp = diffusion_nn(adata_sub, k=k0).astype("float")
                    # call kBET
                    score = kBET_single(
                        matrix=matrix,
                        batch=adata_sub.obs[batch_key],
                        knn=nn_index_tmp
                        + 1,  # nn_index in python is 0-based and 1-based in R
                        verbose=verbose,
                        k0=k0,
                    )
                except NeighborsError:
                    print("Not enough neighbours")
                    score = 1  # i.e. 100% rejection

            else:
                # check the number of components where kBET can be computed upon
                comp_size = pd.value_counts(labs)
                # check which components are small
                comp_size_thresh = 3 * k0
                idx_nonan = np.flatnonzero(
                    np.in1d(labs, comp_size[comp_size >= comp_size_thresh].index)
                )

                # check if 75% of all cells can be used for kBET run
                if len(idx_nonan) / len(labs) >= 0.75:
                    # create another subset of components, assume they are not visited in a diffusion process
                    adata_sub_sub = adata_sub[idx_nonan, :].copy()
                    nn_index_tmp = np.empty(shape=(adata_sub.n_obs, k0))
                    nn_index_tmp[:] = np.nan

                    try:
                        nn_index_tmp[idx_nonan] = diffusion_nn(
                            adata_sub_sub, k=k0
                        ).astype("float")
                        # call kBET
                        score = kBET_single(
                            matrix=matrix,
                            batch=adata_sub.obs[batch_key],
                            knn=nn_index_tmp
                            + 1,  # nn_index in python is 0-based and 1-based in R
                            verbose=verbose,
                            k0=k0,
                        )
                    except NeighborsError:
                        print("Not enough neighbours")
                        score = 1  # i.e. 100% rejection
                else:  # if there are too many too small connected components, set kBET score to 1
                    score = 1  # i.e. 100% rejection

        kBET_scores["cluster"].append(clus)
        kBET_scores["kBET"].append(score)

    kBET_scores = pd.DataFrame.from_dict(kBET_scores)
    kBET_scores = kBET_scores.reset_index(drop=True)

    if return_df:
        return kBET_scores

    final_score = np.nanmean(kBET_scores["kBET"])
    return 1 - final_score if scaled else final_score


def kBET_single(matrix, batch, k0=10, knn=None, verbose=False):
    """Single kBET run

    Compute k-nearest neighbour batch effect test (kBET) score as described in
    https://doi.org/10.1038/s41592-018-0254-1

    :param matrix: expression matrix (at the moment: a PCA matrix, so ``do.pca`` is set to ``FALSE``)
    :param batch: series or list of batch assignments
    :returns: kBET observed rejection rate
    """
    try:
        import anndata2ri
        import rpy2.rinterface_lib.callbacks
        import rpy2.rinterface_lib.embedded
        import rpy2.robjects as ro

        rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    except ModuleNotFoundError as e:
        print("---------")
        raise Exception("rpy2 not available: " + e)


    try:
        ro.r("library(kBET)")
    except rpy2.rinterface_lib.embedded.RRuntimeError as ex:
        raise Exception("Have problem loading package: "+ ex)

    anndata2ri.activate()

    if verbose:
        print("importing expression matrix")
    ro.globalenv["data_mtrx"] = matrix
    ro.globalenv["batch"] = batch

    if verbose:
        print("kBET estimation")

    ro.globalenv["knn_graph"] = knn
    ro.globalenv["k0"] = k0
    ro.r(
        "batch.estimate <- kBET("
        "  data_mtrx,"
        "  batch,"
        "  knn=knn_graph,"
        "  k0=k0,"
        "  plot=FALSE,"
        "  do.pca=FALSE,"
        "  heuristic=FALSE,"
        "  adapt=FALSE,"
        f"  verbose={str(verbose).upper()}"
        ")"
    )

    try:
        score = ro.r("batch.estimate$summary$kBET.observed")[0]
    except rpy2.rinterface_lib.embedded.RRuntimeError as ex:
        print(f"Error computing kBET: {ex}\nSetting value to np.nan")
        score = np.nan

    anndata2ri.deactivate()

    return score