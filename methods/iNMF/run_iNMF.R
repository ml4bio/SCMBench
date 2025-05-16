#!/usr/bin/env Rscript


source("/mnt/nas/user/yixuan/deepmaps-master/scRNA_scATAC1.r")
Sys.setenv(RETICULATE_PYTHON = "/mnt/nas/user/yixuan/miniconda3/envs/R4+torgeo/bin/python")
use_python("/mnt/nas/user/yixuan/miniconda3/envs/R4+torgeo/bin/python")
py_config()
library(argparse)
library(rliger)
library(yaml)
nonneg <- function(x, eps = 1e-16) {
  x[x < eps] = eps
  return(x)
}
# source("/mnt/nas/user/yixuan/Multiomics-benchmark-main/evaluation/workflow/scripts/rliger/rliger.r")
solveNNLS <- function(C, B) {
    .Call('_rliger_solveNNLS', PACKAGE = 'rliger', C, B)
}
online_iNMF_test <- function(object,
                        X_new = NULL,
                        projection = FALSE,
                        W.init = NULL,
                        V.init = NULL,
                        H.init = NULL,
                        A.init = NULL,
                        B.init = NULL,
                        k = 20,
                        lambda = 5,
                        max.epochs = 5,
                        miniBatch_max_iters = 1,
                        miniBatch_size = 5000,
                        h5_chunk_size = 1000,
                        seed = 123,
                        verbose = TRUE){
  if (!is.null(X_new)){ # if there is new dataset
    raw.data_prev = object@raw.data
    norm.data_prev = object@norm.data
    h5file.info_prev = object@h5file.info
    scale.data_prev = object@scale.data
    cell.data_prev = object@cell.data
    names(raw.data_prev) = names(object@raw.data)

    # assuming only one new dataset arrives at a time
    raw.data = c()
    norm.data = c()
    h5file.info = c()
    scale.data = c()
    cell.data = c()
    for (i in 1:length(X_new)){
      raw.data = c(raw.data, X_new[[i]]@raw.data)
      norm.data = c(norm.data, X_new[[i]]@norm.data)
      h5file.info = c(h5file.info, X_new[[i]]@h5file.info)
      scale.data = c(scale.data, X_new[[i]]@scale.data)
      cell.data = rbind(cell.data, X_new[[i]]@cell.data)
    }
    object@raw.data = raw.data
    object@norm.data = norm.data
    object@h5file.info = h5file.info
    object@scale.data = scale.data
    object@cell.data = cell.data

    # check whether X_new needs to be processed
    for (i in 1:length(object@raw.data)){
      if (class(object@raw.data[[i]])[1] == "H5File"){
        processed = object@raw.data[[i]]$exists("scale.data")
      } else {
        processed = !is.null(X_new[[i]]@scale.data)
      }

      if (processed) {
        if (verbose) {
          cat("New dataset", i, "already preprocessed.", "\n")
        }
      } else {
        if (verbose) {
          cat("New dataset", i, "not preprocessed. Preprocessing...", "\n")
        }
        object = normalize(object, chunk = h5_chunk_size)
        object = scaleNotCenter(object, remove.missing = TRUE, chunk = h5_chunk_size)
        if (verbose) {
          cat("New dataset", i, "Processed.", "\n")
        }
      }
    }


    object@raw.data = c(raw.data_prev, object@raw.data)
    object@norm.data = c(norm.data_prev, object@norm.data)
    object@h5file.info = c(h5file.info_prev, object@h5file.info)
    object@scale.data = c(scale.data_prev, object@scale.data)
    object@cell.data = rbind(cell.data_prev, object@cell.data)
    # k x gene -> gene x k & cell x k-> k x cell
    object@W = t(object@W)
    object@V = lapply(object@V, t)
    object@H = lapply(object@H, t)
  }

  for (i in 1:length(object@raw.data)){
    if (class(object@raw.data[[i]])[1] != "H5File") object@scale.data[[i]] = t(object@scale.data[[i]])
  }

  ## extract required information and initialize algorithm
  num_files = length(object@raw.data) # number of total input hdf5 files
  num_prev_files = 0 # number of input hdf5 files processed in last step
  num_new_files = 0 # number of new input hdf5 files since last step
  if (is.null(X_new)) {
    num_prev_files = 0 # start from scratch
    num_new_files = num_files
  } else {
    num_new_files = length(X_new)
    num_prev_files = num_files - num_new_files
    if (verbose) {
      cat(num_new_files, "new datasets detected.", "\n")
    }
  }

  file_idx = 1:num_files # indices for all input files
  file_idx_new = (num_prev_files+1):num_files # indices only for new input files
  file_idx_prev = setdiff(file_idx,file_idx_new)

  vargenes = object@var.genes
  file_names = names(object@raw.data)
  gene_names = vargenes # genes selected for analysis
  num_genes = length(vargenes) # number of the selected genes

  cell_barcodes = list() # cell barcodes for each dataset
  for (i in file_idx){
    cell_barcodes[[i]] = rownames(object@cell.data)[object@cell.data$dataset == file_names[i]]
  }
  num_cells = unlist(lapply(cell_barcodes, length)) # number of cells in each dataset
  num_cells_new = num_cells[(num_prev_files+1):num_files]
  minibatch_sizes = rep(0, num_files)

  for (i in file_idx_new) {
    minibatch_sizes[i] = round((num_cells[i]/sum(num_cells[file_idx_new])) * miniBatch_size)
  }
  minibatch_sizes_orig = minibatch_sizes

  if (!projection) {

    if(!is.null(seed)){
      set.seed(seed)
    }

    # W matrix initialization
    if (is.null(X_new)) {
      object@W = matrix(abs(runif(num_genes * k, 0, 2)), num_genes, k)
      for (j in 1:k){
        object@W[, j] = object@W[, j] / sqrt(sum(object@W[, j]^2))
      }
    } else {
      object@W = if(!is.null(W.init)) W.init else object@W
    }
    # V_i matrix initialization
    if (is.null(X_new)) {
      object@V = list()
      for (i in file_idx){
        V_init_idx = sample(1:num_cells_new[i], k) # pick k sample from datasets as initial H matrix
        # object@V[[i]] = object@scale.data[[i]][1:num_genes, V_init_idx]

        dense_temp = as.matrix(object@scale.data[[i]])
        object@V[[i]] = dense_temp[1:num_genes, V_init_idx]

        #object@V[[i]] = matrix(data = abs(x = runif(n = num_genes * k, min = 0, max = 2)),
        #                       nrow = num_genes,
        #                       ncol = k)
      }
      

      # normalize the columns of H_i, H_s matrices
      for (j in 1:k){
        for (i in file_idx){ # normalize columns of dictionaries
          object@V[[i]][, j] = object@V[[i]][, j] / sqrt(sum(object@V[[i]][, j]^2))
        }
      }
    } else { # if previous Vs are provided
      object@V[file_idx_prev] = if(!is.null(V.init)) V.init else object@V
      V_init_idx = list()
      for (i in file_idx_new){
        V_init_idx = sample(1:num_cells[i], k)
        object@V[[i]] = object@scale.data[[i]][1:num_genes, V_init_idx] # initialize the Vi for new dataset
        for (j in 1:k){
          object@V[[i]][, j] = object@V[[i]][, j] / sqrt(sum(object@V[[i]][, j]^2))
        }
      }
    }
    # H_i matrices initialization
    if (is.null(X_new)) {
      object@H = rep(list(NULL),num_files)
      H_minibatch = list()
    } else { # if previous Hs are provided
      object@H[file_idx_prev] = if(!is.null(H.init)) H.init else object@H
      object@H[file_idx_new] = rep(list(NULL),num_new_files)
      H_minibatch = list()
    }
    # A = HiHi^t, B = XiHit
    A_old = list()
    B_old = list()

    if (is.null(X_new)) {
      object@A = rep(list(matrix(0, k, k)), num_new_files)
      object@B = rep(list(matrix(0, num_genes, k)), num_new_files)
      A_old = rep(list(matrix(0, k, k)), num_new_files) # save information older than 2 epochs
      B_old = rep(list(matrix(0, num_genes, k)), num_new_files) # save information older than 2 epochs

    } else {
      object@A[file_idx_prev] = if(!is.null(A.init)) A.init else object@A
      object@B[file_idx_prev] = if(!is.null(B.init)) B.init else object@B
      A_old[file_idx_prev] = rep(list(NULL), num_prev_files)
      B_old[file_idx_prev] = rep(list(NULL), num_prev_files)
      object@A[(num_prev_files+1):num_files] = rep(list(matrix(0, k, k)), num_new_files)
      object@B[(num_prev_files+1):num_files] = rep(list(matrix(0, num_genes, k)), num_new_files)
      A_old[(num_prev_files+1):num_files] = rep(list(matrix(0, k, k)), num_new_files) # save information older than 2 epochs
      B_old[(num_prev_files+1):num_files] = rep(list(matrix(0, k, k)), num_new_files) # save information older than 2 epochs
    }

    iter = 1
    epoch = rep(0, num_files) # intialize the number of epoch for each dataset
    epoch_prev = rep(0, num_files) # intialize the previous number of epoch for each dataset
    epoch_next = rep(FALSE, num_files)
    sqrt_lambda = sqrt(lambda)
    total_time = 0 # track the total amount of time used for the online learning


    num_chunks = rep(NULL, num_files)
    chunk_idx = rep(list(NULL), num_files)
    all_idx = rep(list(NULL), num_files)

    # chunk permutation
    for (i in file_idx_new){
      num_chunks[i] = ceiling(num_cells[i]/h5_chunk_size)
      chunk_idx[[i]] = sample(1:num_chunks[i],num_chunks[i])
      # idx in the first chunk
      if(chunk_idx[[i]][1]!=num_chunks[i]){
        all_idx[[i]] = (1+h5_chunk_size*(chunk_idx[[i]][1]-1)):(chunk_idx[[i]][1]*h5_chunk_size)
      } else {
        all_idx[[i]] = (1+h5_chunk_size*(chunk_idx[[i]][1]-1)):(num_cells[i])
      }

      for (j in chunk_idx[[i]][-1]){
        if (j != num_chunks[i]){
          all_idx[[i]] = c(all_idx[[i]],(1+h5_chunk_size*(j-1)):(j*h5_chunk_size))
        } else {
          all_idx[[i]] = c(all_idx[[i]],(1+h5_chunk_size*(j-1)):num_cells[i])
        }
      }
    }

    total.iters = floor(sum(num_cells_new) * max.epochs / miniBatch_size)
    if (verbose) {
      cat("Starting Online iNMF...", "\n")
      pb <- txtProgressBar(min = 1, max = total.iters+1, style = 3)
    } 
    
    while(epoch[file_idx_new[1]] < max.epochs) {
      # track epochs
      minibatch_idx = rep(list(NULL), num_files) # indices of samples in each dataest used for this iteration
      if ((max.epochs * num_cells_new[1] - (iter-1) * minibatch_sizes[file_idx_new[1]]) >= minibatch_sizes[file_idx_new[1]]){ # check if the size of the last mini-batch == pre-specified mini-batch size
        for (i in file_idx_new){
          epoch[i] = (iter * minibatch_sizes[i]) %/% num_cells[i] # caculate the current epoch
          if ((epoch_prev[i] != epoch[i]) & ((iter * minibatch_sizes[i]) %% num_cells[i] != 0)){ # if current iter cycles through the data and start a new cycle
            epoch_next[i] = TRUE
            epoch_prev[i] = epoch[i]
            # shuffle dataset before the next epoch
            minibatch_idx[[i]] = all_idx[[i]][c(((((iter - 1) * minibatch_sizes[i]) %% num_cells[i]) + 1):num_cells[i])]
            chunk_idx[[i]] = sample(1:num_chunks[i],num_chunks[i])
            all_idx[[i]] = 0
            for (j in chunk_idx[[i]]){
              if (j != num_chunks[i]){
                all_idx[[i]] = c(all_idx[[i]],(1+h5_chunk_size*(j-1)):(j*h5_chunk_size))
              }else{
                all_idx[[i]] = c(all_idx[[i]],(1+h5_chunk_size*(j-1)):num_cells[i])
              }
            }
            all_idx[[i]] = all_idx[[i]][-1] # remove the first element 0
            minibatch_idx[[i]] = c(minibatch_idx[[i]],all_idx[[i]][1:((iter * minibatch_sizes[i]) %% num_cells[i])])

          } else if ((epoch_prev[i] != epoch[i]) & ((iter * minibatch_sizes[i]) %% num_cells[i] == 0)){ # if current iter finishes this cycle without start a a new cycle
            epoch_next[i] = TRUE
            epoch_prev[i] = epoch[i]

            minibatch_idx[[i]] = all_idx[[i]][((((iter-1) * minibatch_sizes[i]) %% num_cells[i]) + 1):num_cells[i]]
            chunk_idx[[i]] = sample(1:num_chunks[i],num_chunks[i])
            all_idx[[i]] = 0
            for (j in chunk_idx[[i]]){
              if (j != num_chunks[i]){
                all_idx[[i]] = c(all_idx[[i]],(1+h5_chunk_size*(j-1)):(j*h5_chunk_size))
              }else{
                all_idx[[i]] = c(all_idx[[i]],(1+h5_chunk_size*(j-1)):num_cells[i])
              }
            }
            all_idx[[i]] = all_idx[[i]][-1] # remove the first element 0
          } else {                                                                        # if current iter stays within a single cycle
            minibatch_idx[[i]] = all_idx[[i]][(((iter-1) * minibatch_sizes[i]) %% num_cells[i] + 1):((iter * minibatch_sizes[i]) %% num_cells[i])]
          }
        }
      } else {
        for (i in file_idx_new){
          minibatch_sizes[i] = max.epochs * num_cells[i] - (iter-1) * minibatch_sizes[i]
          minibatch_idx[[i]] = (((iter-1) * minibatch_sizes_orig[i] + 1) %% num_cells[i]):num_cells[i]
        }
        epoch[file_idx_new[1]] = max.epochs # last epoch
      }


      if (length(minibatch_idx[[file_idx_new[1]]]) == minibatch_sizes_orig[file_idx_new[1]]){
        X_minibatch = rep(list(NULL), num_files)
        for (i in file_idx_new){
          X_minibatch[[i]] = object@scale.data[[i]][1:num_genes ,minibatch_idx[[i]]]
        }
        # update H_i by ANLS Hi_minibatch[[i]]
        H_minibatch = rep(list(NULL), num_files)
        for (i in file_idx_new){
          H_minibatch[[i]] = solveNNLS(rbind(object@W + object@V[[i]], sqrt_lambda * object@V[[i]]),
                                       as.matrix(rbind(X_minibatch[[i]], matrix(0, num_genes, minibatch_sizes[i]))))
        }
        
        
        # updata A and B matrices
        if (iter == 1){
          scale_param = c(rep(0, num_prev_files), rep(0, num_new_files))
        } else if(iter == 2){
          scale_param = c(rep(0, num_prev_files), rep(1, num_new_files) / minibatch_sizes[file_idx_new])
        } else {
          scale_param = c(rep(0, num_prev_files), rep((iter - 2) / (iter - 1), num_new_files))
        }


        if (epoch[file_idx_new[1]] > 0 & epoch_next[file_idx_new[1]] == TRUE){ # remove information older than 2 epochs
          for (i in file_idx_new){
            object@A[[i]] = object@A[[i]] - A_old[[i]]
            A_old[[i]] = scale_param[i] * object@A[[i]]
            object@B[[i]] = object@B[[i]] - B_old[[i]]
            B_old[[i]] = scale_param[i] * object@B[[i]]
          }
        } else{ # otherwise scale the old information
          for (i in file_idx_new){
            A_old[[i]] = scale_param[i] * A_old[[i]]
            B_old[[i]] = scale_param[i] * B_old[[i]]
          }
        }

        for (i in file_idx_new){
          object@A[[i]] = scale_param[i] * object@A[[i]] + H_minibatch[[i]] %*% t(H_minibatch[[i]]) / minibatch_sizes[i]   # HiHit
          diag(object@A[[i]])[diag(object@A[[i]])==0] = 1e-15
          object@B[[i]] = scale_param[i] * object@B[[i]] + X_minibatch[[i]] %*% t(H_minibatch[[i]]) / minibatch_sizes[i]   # XiHit
        }


        # update W, V_i by HALS
        iter_miniBatch = 1
        delta_miniBatch = Inf
        max_iters_miniBatch = miniBatch_max_iters

        while(iter_miniBatch <= max_iters_miniBatch){
          # update W
          for (j in 1:k){
            W_update_numerator = rep(0, num_genes)
            W_update_denominator = 0
            for (i in file_idx){
              W_update_numerator = W_update_numerator + object@B[[i]][, j] - (object@W + object@V[[i]]) %*% object@A[[i]][, j]
              W_update_denominator = W_update_denominator +  object@A[[i]][j,j]
            }

            object@W[, j] = nonneg(object@W[, j] + W_update_numerator / W_update_denominator)
          }

          # update V_i
          for (j in 1:k){
            for (i in file_idx_new){
              object@V[[i]][, j] = nonneg(object@V[[i]][, j] + (object@B[[i]][, j] - (object@W + (1 + lambda) * object@V[[i]]) %*% object@A[[i]][, j]) /
                ((1 + lambda) * object@A[[i]][j, j]))
            }
          }

          iter_miniBatch = iter_miniBatch + 1
        }
        epoch_next = rep(FALSE, num_files) # reset epoch change indicator
        iter = iter + 1
        if (verbose) {
          setTxtProgressBar(pb = pb, value = iter)
        }
      }
    }
    if (verbose) {
      cat("\nCalculate metagene loadings...", "\n")
    }
    object@H = rep(list(NULL), num_files)
    for (i in file_idx){
      if (num_cells[i] %% miniBatch_size == 0) num_batch = num_cells[i] %/% miniBatch_size else num_batch = num_cells[i] %/% miniBatch_size + 1
      if (num_batch == 1){
        X_i = object@scale.data[[i]][1:num_genes,]
        object@H[[i]] = solveNNLS(rbind(object@W + object@V[[i]],sqrt_lambda * object@V[[i]]), as.matrix(rbind(X_i, matrix(0, num_genes , num_cells[i]))))
      } else {
        for (batch_idx in 1:num_batch){
          if (batch_idx != num_batch){
            cell_idx = ((batch_idx - 1) * miniBatch_size + 1):(batch_idx * miniBatch_size)
          } else {
            cell_idx = ((batch_idx - 1) * miniBatch_size + 1):num_cells[i]
          }
          X_i_batch = object@scale.data[[i]][1:num_genes,cell_idx]
          object@H[[i]] = cbind(object@H[[i]], solveNNLS(rbind(object@W + object@V[[i]], sqrt_lambda * object@V[[i]]),
                                                         as.matrix(rbind(X_i_batch, matrix(0, num_genes , length(cell_idx))))))
        }
      }
      colnames(object@H[[i]]) = cell_barcodes[[i]]
    }

    rownames(object@W) = gene_names
    colnames(object@W) = NULL

    for (i in file_idx){
      rownames(object@V[[i]]) = gene_names
      colnames(object@V[[i]]) = NULL
    }

  } else {
    if (verbose) {
      cat("Metagene projection", "\n")
    }
    object@W = if(!is.null(W.init)) W.init else object@W
    object@H[file_idx_new] = rep(list(NULL), num_new_files)
    object@V[file_idx_new] = rep(list(NULL), num_new_files)
    for (i in file_idx_new){
      if (num_cells[i] %% miniBatch_size == 0) num_batch = num_cells[i] %/% miniBatch_size else num_batch = num_cells[i] %/% miniBatch_size + 1
      if (num_cells[i] <= miniBatch_size){
        object@H[[i]] = solveNNLS(object@W, object@scale.data[[i]][1:num_genes,])
      } else {
        for (batch_idx in 1:num_batch){
          if (batch_idx != num_batch){
            cell_idx = ((batch_idx - 1) * miniBatch_size + 1):(batch_idx * miniBatch_size)
          } else {
            cell_idx = ((batch_idx - 1) * miniBatch_size + 1):num_cells[i]
          }
          object@H[[i]] = cbind(object@H[[i]],solveNNLS(object@W, object@scale.data[[i]][1:num_genes,cell_idx]))
        }
      }
      colnames(object@H[[i]]) = cell_barcodes[[i]]
      object@V[[i]] = matrix(0, num_genes, k)
    }
  }

  # gene x k -> k x gene & k x cell -> cell x k
  object@W = t(object@W)
  object@V = lapply(object@V, t)
  object@H = lapply(object@H, t)
  for (i in 1:length(object@raw.data)){
    if (class(object@raw.data[[i]])[1] != "H5File") object@scale.data[[i]] = t(object@scale.data[[i]])
  }

  if (!is.null(X_new)){
    names(object@scale.data) <- names(object@raw.data) <- c(names(raw.data_prev), names(X_new))
  }
  names(object@H) <- names(object@V) <- names(object@raw.data)
  return(object)
}

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
    # hvg <- intersect(hvg, rownames(int.liger@raw.data$rna))  # Because of feature filtering in createLiger
    # hvg <- intersect(hvg, rownames(int.liger@raw.data$atac))  # Because of feature filtering in createLiger
    rna_cells <- rownames(rna$obs)
    atac_cells <- rownames(atac$obs)
    n_cells <- nrow(rna$obs) + nrow(atac$obs)
    min_cells <- min(nrow(rna$obs), nrow(atac$obs))
    # print(rna, atac)
    rm(rna, atac)
    gc()  # Reduce memory usage
    cat("[2/4] Data preprocessing...\n")
    start_time <- proc.time()
    int.liger <- normalize(int.liger,remove.missing = FALSE)
    int.liger <- scaleNotCenter(int.liger,remove.missing = FALSE)
    int.liger@scale.data$rna <- t(int.liger@norm.data$rna)
    int.liger@scale.data$atac <- t(int.liger@norm.data$atac)
    int.liger@var.genes <- hvg
    cat("[3/4] Running iNMF...\n")
    int.liger <- online_iNMF_test(int.liger, k = 20, miniBatch_size = min(5000, min_cells), seed = args$random_seed)
    int.liger <- quantile_norm(int.liger, rand.seed = args$random_seed)
    elapsed_time <- proc.time() - start_time

    cat("[4/4] Saving results...\n")
    combined_latent <- int.liger@H.norm
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
    atac_latent <- combined_latent[atac_cells, ]
    rownames(atac_latent) <- gsub("\\.ATAC$", "", rownames(atac_latent))
    write.table(
        rna_latent, args$output_rna,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
    )
    write.table(
        atac_latent, args$output_atac,
        sep = ",", row.names = TRUE, col.names = FALSE, quote = FALSE, na = ""
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
