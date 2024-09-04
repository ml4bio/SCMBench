library(Seurat)
options(Seurat.object.assay.version = "v5")
library(Signac)
library(EnsDb.Hsapiens.v86)
library(cowplot)
library(Matrix)
library(dplyr)
library(yaml)
library(argparse)
library(ggplot2)

library(stringr)




parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument(
        "--input-bridge-rna", dest = "input_bridge_rna", type = "character", required = TRUE,
        help = "Path to input bridge RNA dataset (.h5ad)"
    )
    parser$add_argument(
        "--input-bridge-atac", dest = "input_bridge_atac", type = "character", required = TRUE,
        help = "Path to input bridge ATAC dataset (.h5ad)"
    )
    parser$add_argument(
        "--input-bridge-fragments", dest = "input_bridge_fragments", type = "character", required = TRUE,
        help = "Path to input fragments (.tsv.gz)"
    )
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

    rna_bridge<-read_h5ad(args$input_bridge_rna)$X
    atac_bridge<-read_h5ad(args$input_bridge_atac)$X

    stopifnot(dim(atac_bridge)[1] == dim(rna_bridge)[1])

    rna<-read_h5ad(args$input_rna)$X
    atac<-read_h5ad(args$input_atac)$X

    cat("[2/4] Data preprocessing...\n")
    start_time <- proc.time()
    
    # construct bridge object
    obj.multi <- CreateSeuratObject(counts = t(rna_bridge))
    obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")
    
    #grange.counts <- StringToGRanges(rownames(atac_bridge), sep = c(":", "-"))
    #grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
    seqlevels(annotations) <- ucsc.levels
    genome(annotations) <- "hg38"

    bridge_frag.file <- args$input_bridge_fragments
    chrom_assay <- CreateChromatinAssay(counts = t(atac_bridge), sep = c(":", "-"), genome = 'hg38', 
    fragments = bridge_frag.file,
    min.cells = 10,
    annotation = annotations)

    obj.multi[["ATAC"]] <- chrom_assay

    obj.multi <- subset(
    x = obj.multi,
    subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
    )
    # construct ATAC object

    fragpath <- args$input_fragments


    atac_assay <- CreateChromatinAssay(
    counts = t(atac),
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotations
    )

    requant_multiome_ATAC <- FeatureMatrix(
    fragments = Fragments(atac_assay),
    features = granges(obj.multi[['ATAC']]),
    cells = Cells(atac_assay)
    )

    ATAC_assay <- CreateChromatinAssay(
    counts = requant_multiome_ATAC,
    fragments = fragpath,
    annotation = annotations
    )

    obj.atac  <- CreateSeuratObject(counts = ATAC_assay,assay = 'ATAC')
    obj.atac[['peak.orig']] <- atac_assay
    obj.atac <- subset(obj.atac, subset = nCount_ATAC < 7e4 & nCount_ATAC > 2000)
    
    obj.rna <- CreateSeuratObject(counts=t(rna))
    obj.rna = SCTransform(object = obj.rna) %>%
              RunPCA() %>%
              RunUMAP(dims = 1:50, return.model = TRUE)
    #Preprocessing/normalization for all datasets
    # normalize multiome RNA
    DefaultAssay(obj.multi) <- "RNA"
    obj.multi <- SCTransform(obj.multi, verbose = FALSE)
    # normalize multiome ATAC
    DefaultAssay(obj.multi) <- "ATAC"
    obj.multi <- RunTFIDF(obj.multi)
    obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
    # normalize query
    obj.atac <- RunTFIDF(obj.atac)

    cat("[3/4] Running Seurat v5...\n")
    # Drop first dimension for ATAC reduction
    dims.atac <- 2:50
    dims.rna <- 1:50
    DefaultAssay(obj.multi) <-  "RNA"
    DefaultAssay(obj.rna) <- "SCT"
    obj.rna.ext <- PrepareBridgeReference(
    reference = obj.rna, bridge = obj.multi,
    reference.reduction = "pca", reference.dims = dims.rna,
    normalization.method = "SCT")

    bridge.anchor <- FindBridgeTransferAnchors(
    extended.reference = obj.rna.ext, query = obj.atac,
    reduction = "lsiproject", dims = dims.atac)

    obj.atac <- MapQuery(
    anchorset = bridge.anchor, reference = obj.rna.ext,
    query = obj.atac,
    reduction.model = "umap")

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells

    rna_latent<-obj.rna.ext@reductions$pca@cell.embeddings

    atac_latent<-obj.atac@reductions$ref.Bridge.reduc@cell.embeddings


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
    write_yaml(
        list(
            args = args,
            time = elapsed_time["elapsed"],
            rna_cells = dim(obj.rna)[2],
            atac_cells = dim(obj.atac)[2]
        ), args$run_info
    )
}

main(parse_args())
