#!/usr/bin/env Rscript

source("../../SCMBench/scATAC_source.r")
Sys.setenv(RETICULATE_PYTHON = "/mnt/nas/user/yixuan/miniconda3/envs/R4+torgeo/bin/python")
use_python("/mnt/nas/user/yixuan/miniconda3/envs/R4+torgeo/bin/python")
py_config()
library(argparse)
library(Seurat)
library(yaml)


parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument(
        "--input-rna", dest = "input_rna", type = "character", required = TRUE,
        help = "Path to input RNA dataset (.h5ad)"
    )
    parser$add_argument(
        "--input-atac", dest = "input_atac", type = "character", required = TRUE,
        help = "Path to input ATAC dataset converted to RNA space (.h5ad)"
    )
    parser$add_argument(
        "-s", "--random-seed", dest = "random_seed", type = "integer", default = 0,
        help = "Random seed"
    )
    parser$add_argument(
        "--output-active", dest = "output_active", type = "character", required = TRUE,
        help = "Path of output ATAC activity file (.h5ad)"
    )
    parser$add_argument(
        "--genome", dest = "genome", type = "character", required = TRUE,
        help = "mouse or human"
    )
    return(parser$parse_args())
}

read_h5ad <- function(filename) {
    builtins <- reticulate::import_builtins()
    anndata <- reticulate::import("anndata", convert = FALSE)

    Mapping <- reticulate::import("typing")$Mapping
    DataFrame <- reticulate::import("pandas")$DataFrame
    issparse <- reticulate::import("scipy.sparse")$issparse
    isinstance <- builtins$isinstance

    adata <- anndata$read_h5ad(filename)

    .convert <- function(obj) {
        if (!isinstance(obj, Mapping) || isinstance(obj, DataFrame)) {
            return(reticulate::py_to_r(obj))
        }
        ret <- list()
        for (item in builtins$list(obj$keys())) {
            ret[[item]] <- .convert(obj[[item]])
        }
        return(ret)
    }

    if (issparse(adata$X)) {
        X <- .convert(adata$X$tocsc())
    } else {
        X <- .convert(adata$X)
    }
    layers <- .convert(adata$layers)
    obs <- .convert(adata$obs)
    var <- .convert(adata$var)
    obsm <- .convert(adata$obsm)
    varm <- .convert(adata$varm)
    obsp <- .convert(adata$obsp)
    varp <- .convert(adata$varp)
    uns <- .convert(adata$uns)
    rownames(X) <- rownames(obs)
    colnames(X) <- rownames(var)

    return(list(
        X = X, layers = layers,
        obs = obs, var = var,
        obsm = obsm, varm = varm,
        obsp = obsp, varp = varp,
        uns = uns
    ))
}


main <- function(args) {

    set.seed(args$random_seed)

    cat("[1/3] Reading data...\n")
    rna <- read_h5ad(args$input_atac)
    atac <- read_h5ad(args$input_atac)
    # stopifnot(all(rownames(atac$var) == rownames(rna$var)))
    rownames(rna$obs) <- paste(rownames(rna$obs))  # Avoid collision
    rownames(rna$X) <- rownames(rna$obs)
    colnames(rna$X) <- rownames(rna$var)
    rownames(atac$obs) <- paste(rownames(atac$obs))  # Avoid collision
    rownames(atac$X) <- rownames(atac$obs)
    colnames(atac$X) <- rownames(atac$var)
    print(dim(t(rna$X)))
    print(dim(t(atac$X)))
    lymph_obj <- ReadData(h5Path = Null, rna_matrix = t(rna$X), atac_matrix =t(atac$X), data_type = "scRNA_scATAC", dataFormat = "matrixs", cell.filter=FALSE, min_cell = 0.1)
    # print(dim(lymph_obj))
    print(class(lymph_obj@assays$ATAC@counts))
    print(class(atac$X))
    # rna.so <- CreateSeuratObject(counts = Matrix::t(rna$X), meta.data = rna$obs, min_cell = 0.001)
    # atac.so <- CreateSeuratObject(counts = Matrix::t(atac$X), meta.data = atac$obs, min_cell = 0.001)

    hvg <- rownames(rna$var)[rna$var$highly_variable]
    n_cells <- nrow(rna$obs) + nrow(atac$obs)


    cat("[2/3] Generating Activity...\n")
    activity_matrix <- ATACCalculateGenescore(lymph_obj@assays$ATAC@counts,organism = args$genome)
    # activity_matrix <- ATACCalculateGenescore(atac$X,organism = args$genome)
    activity<-t(activity_matrix)
    library(anndata)
    activity_ad <- AnnData(X = activity,obs = atac$obs,var = colnames(activity))
    activity_ad$var<-data.frame(colnames(activity), row.names = colnames(activity))
    
    cat("[3/3] Saving results...\n")
    write_h5ad(activity_ad, args$output_active)

}


main(parse_args())
