#!/usr/bin/env Rscript

library(reticulate)
use_condaenv(condaenv = "/workspace/wangyixuan/.conda/envs/zgy_py3.8", required = TRUE)
library(jsonlite)            
comp <- import_from_path("comp_effic", path = "/mnt/nas/user/yixuan/zgy/SCMBench")

library(argparse)
library(rliger)
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
        "--device", dest="device", type="character", default='cpu',
        help="cuda or cpu"
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
    rownames(rna$obs) <- paste(rownames(rna$obs), "RNA", sep = ".")  # Avoid collision
    rownames(rna$X) <- rownames(rna$obs)
    colnames(rna$X) <- rownames(rna$var)
    rownames(atac$obs) <- paste(rownames(atac$obs), "ATAC", sep = ".")  # Avoid collision
    rownames(atac$X) <- rownames(atac$obs)
    colnames(atac$X) <- rownames(atac$var)
    # print(head(rna$obs))
    # print(head(atac$obs))
    # print(head(rna$var))
    # print(head(atac$var))
    # print(dim(rna$X))
    # print(dim(atac$X))
    if (args$device == 'cpu') {
        run_info <- comp$log_system_info(use_gpu = FALSE)
        run_info$r_version <- R.version$version.string
        monitor <- comp$PerformanceMonitor("iNMF_cpu", use_gpu = FALSE, gpu_device = 0)
    } else {
        run_info <- comp$log_system_info(use_gpu = TRUE)
        run_info$r_version <- R.version$version.string
        monitor <- comp$PerformanceMonitor("iNMF_gpu", use_gpu = TRUE, gpu_device = 0)
    }
    monitor$start()

    int.liger <- createLiger(list(
        atac = Matrix::t(atac$X),
        rna = Matrix::t(rna$X)
    ))
    # print(int.liger@scale.data)
    # int.liger <- createLiger(list(
    #     atac = atac$X,
    #     rna = rna$X
    # ))
    # print(head(rna$var))
    # print(head(atac$var))
    hvg <- rownames(rna$var)[rna$var$highly_variable]
    # stopifnot(all(hvg %in% rownames(rna$var)))
    # stopifnot(all(hvg %in% rownames(atac$var)))
    hvg <- intersect(hvg, rownames(int.liger@datasets$rna))  # Because of feature filtering in createLiger
    hvg <- intersect(hvg, rownames(int.liger@datasets$atac))  # Because of feature filtering in createLiger
    rna_cells <- rownames(rna$obs)
    atac_cells <- rownames(atac$obs)
    n_cells <- nrow(rna$obs) + nrow(atac$obs)
    min_cells <- min(nrow(rna$obs), nrow(atac$obs))
    # print(rna, atac)
    rm(rna, atac)
    gc()  # Reduce memory usage
    cat("[2/4] Data preprocessing...\n")


    # start_time <- proc.time()
    int.liger <- normalize(int.liger,remove.missing = FALSE)
    int.liger@varFeatures <- hvg
    # int.liger <- selectGenes(int.liger, useDatasets = "rna")  
    int.liger <- scaleNotCenter(int.liger,remove.missing = FALSE)

    cat("[3/4] Running iNMF...\n")
    int.liger <- online_iNMF(int.liger, k = 20, miniBatch_size = min(5000, min_cells), seed = args$random_seed) #### iNMF

    int.liger <- quantile_norm(int.liger, rand.seed = args$random_seed)

    cat("[4/4] Saving results...\n")
    combined_latent <- int.liger@H.norm
    rna_cells = paste0('rna_',rna_cells)
    atac_cells = paste0('atac_',atac_cells)
    missing_cells <- setdiff(
        union(rna_cells, atac_cells),
        rownames(combined_latent)
    )  # Because of cell filtering in scaleNotCenter
    combined_latent <- rbind(combined_latent, matrix(
        nrow = length(missing_cells),
        ncol = ncol(combined_latent),
        dimnames = list(missing_cells, colnames(combined_latent))
    ))  # Fill with NA
    rna_latent <- combined_latent[rna_cells, ]
    rownames(rna_latent) <- gsub("\\.RNA$", "", rownames(rna_latent))
    rownames(rna_latent) <- gsub("rna_", "", rownames(rna_latent))

    atac_latent <- combined_latent[atac_cells, ]
    rownames(atac_latent) <- gsub("\\.ATAC$", "", rownames(atac_latent))
    rownames(atac_latent) <- gsub("atac_", "", rownames(atac_latent))
    # elapsed_time <- proc.time() - start_time
    stats <- monitor$stop()
    run_info <- c(run_info, stats)
    write_json(run_info, path = args$run_info, pretty = TRUE, auto_unbox = TRUE)

    write.table(
        rna_latent, args$output_rna,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
    )
    write.table(
        atac_latent, args$output_atac,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
    )
    # write_yaml(
    #     list(
    #         args = args,
    #         time = elapsed_time["elapsed"],
    #         n_cells = n_cells
    #     ), args$run_info
    # )
}


main(parse_args())
