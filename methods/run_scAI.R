library(argparse)
library(scAI)
library(yaml)
library(httr)
library(dplyr)
library(anndata)
library(reticulate)
library(Matrix)

httr::set_config(httr::config(ssl_verifypeer = FALSE))

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
        "--output", dest = "output_atac", type = "character", required = TRUE,
        help = "Path of output ATAC latent file (.csv)"
    )
    parser$add_argument(
        "--run-info", dest = "run_info", type = "character", required = TRUE,
        help = "Path of output run info file (.yaml)"
    )
    return(parser$parse_args())
}


main <- function(args) {

    set.seed(args$random_seed)
    cat("[1/4] Reading data...\n")

    anndata <- reticulate::import("anndata", convert = FALSE)
    numpy <- reticulate::import("numpy", convert = FALSE)
    
    X_rna <- anndata$read_h5ad(args$input_rna)$to_df()
    X_atac <- anndata$read_h5ad(args$input_atac)$to_df()
    labels <- anndata$read_h5ad(args$input_rna)$obs
    X_rna <- py_to_r(X_rna)
    X_rna <- Matrix(as.matrix(X_rna), sparse = TRUE)
    X_rna <- t(X_rna)
    X_atac <- py_to_r(X_atac)
    X_atac <- Matrix(as.matrix(X_atac), sparse = TRUE)
    X_atac <- t(X_atac)
    X <- list(RNA=X_rna, ATAC=X_atac)
    n_cells <- nrow(labels)
    cat("[2/4] Data preprocessing...\n")
    start_time <- proc.time()
    scAI_outs <- create_scAIobject(raw.data = X)

    scAI_outs <- preprocessing(scAI_outs, assay = NULL, minFeatures = 200, minCells = 1, 
                                libararyflag = F, logNormalize = F)

    scAI_outs <- selectFeatures(scAI_outs, assay = "RNA")

    cat("[3/4] Running scAI...\n")
    scAI_outs <- run_scAI(scAI_outs, K = 5, nrun = 5,do.fast = TRUE,hvg.use1 = TRUE)
    elapsed_time <- proc.time() - start_time

    cat("[4/4] Saving results...\n")
    write.table(
        scAI_outs@agg.data, args$output,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE
    )

    write_yaml(
        list(
            args = args,
            time = elapsed_time["elapsed"],
            n_cells = n_cells
        ), args$run_info
    )
    #TBD
}


main(parse_args())
