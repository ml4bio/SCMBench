#!/bin/bash

conda create -n env4 r-base=4.3.3 -y
conda activate env4
conda install conda-forge::r-argparse -y
conda install conda-forge::r-jsonlite -y
conda install conda-forge::r-httr -y
conda install conda-forge::anndata -y
conda install conda-forge::r-reticulate -y
Rscript env4_R_install.r