#!/usr/bin/env Rscript

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
        "--output-rna", dest = "output_rna", type = "character", required = TRUE,
        help = "Path of output RNA latent file (.csv)"
    )
    parser$add_argument(
        "--output-atac", dest = "output_atac", type = "character", required = TRUE,
        help = "Path of output ATAC latent file (.csv)"
    )
    parser$add_argument(
        "--run-info", dest = "run_info", type = "character", required = TRUE,
        help = "Path of output run info file (.yaml)"
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

    cat("[1/4] Reading data...\n")
    rna <- read_h5ad(args$input_rna)
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
    print(dim(lymph_obj))
    
    # rna.so <- CreateSeuratObject(counts = Matrix::t(rna$X), meta.data = rna$obs, min_cell = 0.001)
    # atac.so <- CreateSeuratObject(counts = Matrix::t(atac$X), meta.data = atac$obs, min_cell = 0.001)

    hvg <- rownames(rna$var)[rna$var$highly_variable]
    n_cells <- nrow(rna$obs) + nrow(atac$obs)

    # combined.so <- merge(rna.so, atac.so)
    # rm(rna, atac)
    # gc()  # Reduce memory usage

    cat("[2/4] Data preprocessing...\n")
    start_time <- proc.time()
    # combined.so <- NormalizeData(combined.so)
    # VariableFeatures(combined.so) <- hvg
    # combined.so <- ScaleData(combined.so)
    # combined.so <- RunPCA(combined.so, seed.use = args$random_seed, verbose = FALSE)

    cat("[3/4] Running Deepmaps...\n")
    # activity_matrix <- ATACCalculateGenescore(lymph_obj@assays$ATAC@counts,organism = "GRCm38")
    # activity<-t(activity_matrix)
    # activity_ad <- AnnData(X = activity,obs = atac$obs,var = colnames(activity))
    # activity_ad$var<-data.frame(colnames(activity), row.names = colnames(activity))
    # write_h5ad(activity_ad, "/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/Chen-2019/Chen-2019-activity-ATAC.h5ad")
    # # print(dim(ATAC_gene_peak))
    ATAC_gene_peak <- ATACCalculateGenescore(lymph_obj@assays$ATAC@counts,organism = args$genome)
    # ATAC_gene_peak <- CalGenePeakScore(peak_count_matrix = lymph_obj@assays$ATAC@counts,organism = args$genome)
    GAS_obj <- calculate_GAS_v1(ATAC_gene_peak = ATAC_gene_peak, obj = lymph_obj, method = "wnn", veloPath = NULL)
    GAS <- GAS_obj[[1]]
    lymph_obj <- GAS_obj[[2]]
    HGT_result <- run_HGT(GAS = as.matrix(GAS),result_dir='/mnt/nas/user/yixuan/deepmaps-master/RNA_ATAC', data_type='scRNA_scATAC', envPath=NULL, lr=0.2, epoch=30, n_hid=128, n_heads=16)
    cell_hgt_matrix <- HGT_result[['cell_hgt_matrix']]
    rownames(cell_hgt_matrix) <- colnames(GAS)

    lymph_obj <- lymph_obj[, colnames(GAS)]
    cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]

    # combined.so <- RunHarmony(combined.so, group.by.vars = "domain")
    elapsed_time <- proc.time() - start_time

    cat("[4/4] Saving results...\n")
    # combined_latent <- Embeddings(combined.so, reduction="harmony")
    rna_latent_obj <- CreateDimReducObject(embeddings = cell_hgt_matrix,
                        key = "HGT_",
                        assay = "RNA")
    rna_latent <- Embeddings(rna_latent_obj)
    rownames(rna_latent) <- gsub("\\.RNA$", "", rownames(rna_latent))
    atac_latent_obj <- CreateDimReducObject(embeddings = cell_hgt_matrix,
                        key = "HGT_",
                        assay = "ATAC")
    atac_latent <- Embeddings(atac_latent_obj)
    rownames(atac_latent) <- gsub("\\.ATAC$", "", rownames(atac_latent))
    write.table(
        rna_latent, args$output_rna,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
    )
    write.table(
        atac_latent, args$output_atac,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
    )
    write_yaml(
        list(
            args = args,
            time = elapsed_time["elapsed"],
            n_cells = n_cells
        ), args$run_info
    )
}


main(parse_args())
