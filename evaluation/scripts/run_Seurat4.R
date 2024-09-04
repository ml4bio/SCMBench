library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(yaml)
library(argparse)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(stringr)

parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument(
        "--input-rna", dest = "input_rna", type = "character", required = TRUE,
        help = "Path to input RNA dataset (.h5ad)"
    )
    parser$add_argument(
        "--input-atac", dest = "input_atac", type = "character", required = TRUE,
        help = "Path to input ATAC dataset (.h5ad)"
    )
    parser$add_argument(
        "--input-fragments", dest = "input_fragments", type = "character", required = TRUE,
        help = "Path to input fragments (.tsv.gz)"
    )
    parser$add_argument(
        "-s", "--random-seed", dest = "random_seed", type = "integer", default = 0,
        help = "Random seed"
    )
    parser$add_argument(
        "--output-rna", dest = "output_rna", type = "character", required = TRUE,
        help = "Path of output rna latent file (.csv)"
    )
    parser$add_argument(
        "--output-atac", dest = "output_atac", type = "character", required = TRUE,
        help = "Path of output atac latent file (.csv)"
    )
    parser$add_argument(
        "--output-activity", dest = "output_activity", type = "character", required = TRUE,
        help = "Path of output rna latent file (.csv)"
    )
    parser$add_argument(
        "--run-info", dest = "run_info", type = "character", required = TRUE,
        help = "Path of output run info file (.yaml)"
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

main <- function(args){
    set.seed(args$random_seed)
    cat("[1/4] Reading data...\n")

    pbmc.rna<-read_h5ad(args$input_rna)$X
    pbmc.atac<-read_h5ad(args$input_atac)$X

    cat("[2/4] Data preprocessing...\n")
    start_time <- proc.time()
    pbmc.rna<-CreateSeuratObject(counts=t(pbmc.rna))
    pbmc.atac<-CreateChromatinAssay(counts=t(pbmc.atac),sep = c(":", "-"),genome = 'hg38',fragments=args$input_fragments)
    pbmc.atac<-CreateSeuratObject(counts=pbmc.atac,assay='ATAC')

    pbmc.rna <- NormalizeData(pbmc.rna)
    pbmc.rna <- FindVariableFeatures(pbmc.rna)
    pbmc.rna <- ScaleData(pbmc.rna)
    pbmc.rna <- RunPCA(pbmc.rna)
    pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

    # ATAC analysis add gene annotation information
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    #seqlevelsStyle(annotations) <- "UCSC" some issues in connection
    # substitute steps
    ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
    seqlevels(annotations) <- ucsc.levels
    genome(annotations) <- "hg38"
    Annotation(pbmc.atac) <- annotations

    pbmc.atac <- RunTFIDF(pbmc.atac)
    pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
    pbmc.atac <- RunSVD(pbmc.atac)
    pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

    gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))
    pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
    DefaultAssay(pbmc.atac) <- "ACTIVITY"
    pbmc.atac <- NormalizeData(pbmc.atac)
    pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

    cat("[3/4] Running Seurat v4...\n")
    transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
    genes.use <- VariableFeatures(pbmc.rna)                             
    refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],dims = 2:30)
    pbmc.atac[["RNA"]] <- imputation

    rna_matrix <- GetAssayData(pbmc.rna, assay = "RNA")
    rna_latent<-rna_matrix[genes.use,]

    atac_latent<-GetAssayData(pbmc.atac,assay="RNA")

    activity<-GetAssayData(pbmc.atac,assay="ACTIVITY")

    elapsed_time <- proc.time() - start_time
    cat("[4/4] Saving results...\n")
    write.table(
        t(rna_latent), args$output_rna,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
    )
    write.table(
        t(atac_latent), args$output_atac,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
    )
    write.table(
        t(activity), args$output_activity,
        sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE
    )
    write_yaml(
        list(
            args = args,
            time = elapsed_time["elapsed"],
            n_cells_rna = dim(rna_latent)[2],
            n_cells_atac = dim(atac_latent)[2]
        ), args$run_info
    )
}

main(parse_args())
