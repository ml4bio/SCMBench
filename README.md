# Dev Code for Multiomics Benchmark

Code adapted from [GLUE](https://github.com/gao-lab/GLUE). Major refactoring needed!

Some usage details can also be found in their [docs](https://scglue.readthedocs.io)



## Directory structure

```
.
├── SCMBench                # Main Python package
├── data                    # Data files
├── evaluation              # Method evaluation pipelines
├── methods                 # Packages of several tools
├── run                     # Scripts for running tools
├── env.yaml                # Reproducible Python environment via conda
├── pyproject.toml          # Python package metadata
├── LICENSE
└── README.md
```

Can ignore the `snakemake` for now.

## Quick start
1. Data Preprocessing
Follow [this](https://SCMBench.readthedocs.io/en/latest/preprocessing.html).

2. Run python scrips. E.g., scGLUE.\:
```bash
cd run
mkdir glue-output
python run_GLUE.py --input-rna Chen-2019-RNA.h5ad --input-atac  Chen-2019-ATAC-preprocessed.h5ad -p guidance.graphml.gz --train-dir ./glue-output --output-rna ./glue-output/rna.csv --output-atac ./glue-output/atac.csv --output-feature ./glue-output/features.csv -r glue-output/run_info.yaml
```

3. Run R script. E.g., Harmony.\
Tested approach: install R and related packages with in conda (e.g. [this](https://stackoverflow.com/questions/70410968/is-it-possible-to-install-r-in-miniconda)).
```bash
cd run
mkdir harmony-output
Rscript run_Harmony.R --input-rna Chen-2019-RNA.h5ad --input-atac  Chen-2019-ATAC-preprocessed.h5ad --output-rna ./harmony-output/rna.csv --output-atac ./harmony-output/atac.csv --run-info harmony-output/run_info.yaml
```

## Data preprocessing

Preprocess all the data to the same `.h5ad` format with the same keys following [here](data/README.md).

For new datasets we include, upload preprocessing scripts and update the preprocessed data/link [here](data/README.md) (store the preprocessed data somewhere we can directly download).

Current included datasets:
- ├─ 10x-ATAC-Brain5k\
  ├─ 10x-Multiome-Pbmc10k\
  ├─ Cao-2020\
  ├─ Chen-2019\
  ├─ Domcke-2020\
  ├─ Luo-2017\
  ├─ Ma-2020\
  └─ Saunders-2018
- eqtl\
  └─ GTEx-v8
- hic\
  └─ Javierre-2016
- chip\
  └─ ENCODE\
     └─ TF-human
- database\
   └─ TRRUST-v2

## Algorithms

For newly added downstream tasks, follow the code fashion of `run_[METHOD].py` in [run/scripts/](evaluation/workflow/scripts) and [run/shell/](run/shell/).

Note: Most tools can be directly used by installing it individually via `pip install` if they provide that option. But for scJoint, scGPT, and UCE, we provide edited package in [methods/](methods/) used for their corresponding scripts `run_[METHOD].py`.

Current included algorithms:

Statistical-Based:

Paired:
- [bindSC](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02679-x)
- [MMD_MA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8496402/)
- [MOFA](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)
- [UnionCom](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i48/5870490)

Unpaired:
- [Seurat v5](https://www.nature.com/articles/s41587-023-01767-y)
- [Seurat V4](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3)
- [iNMF](https://www.nature.com/articles/s41587-021-00867-x)
- [LIGER](https://www.cell.com/cell/pdf/S0092-8674(19)30504-5.pdf)
- [Pamona](https://academic.oup.com/bioinformatics/article/38/1/211/6353029)
- [scMoMaT](https://www.nature.com/articles/s41467-023-36066-2)

Deep Learning-based:

Paired:
- [totalVI](https://www.nature.com/articles/s41592-020-01050-x)
- [scMDC](https://www.nature.com/articles/s41467-022-35031-9)
- [DeepMAPS](https://www.nature.com/articles/s41467-023-36559-0)
- [SIMBA](https://www.nature.com/articles/s41592-023-01899-8)

Unpaired:
- [Cobolt](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02556-z)
- [scVI](https://www.nature.com/articles/s41592-018-0229-2)
- [scJoint](https://www.nature.com/articles/s41587-021-01161-6)
- [GLUE](https://www.nature.com/articles/s41587-022-01284-4)
- [Harmony](https://www.nature.com/articles/s41592-019-0619-0)

Foundation Models:
- [Geneformer](https://www.nature.com/articles/s41586-023-06139-9)
- [scFoundation](https://www.nature.com/articles/s41592-024-02305-7)
- [scGPT](https://www.nature.com/articles/s41592-024-02201-0)
- [UCE](https://www.biorxiv.org/content/10.1101/2023.11.28.568918v1)

## Downstream tasks

For newly added downstream tasks, follow the code fashion in [Experiments/](experiments).

Currently included downstream tasks:
- Multi-Omic Integration Accuracy (MAP, MNI, ASW, ARI)
- Bio-Conservation: 
  - Biomarker Detection (JSI)
  - Trajectory Conservation
- Batch effect & detecting over-correction (iLISI, kBET, Graph connectivity, Batch-ASW)

## Solution to potential install issues

- Bedtools. `NotImplementedError: "intersectBed" does not appear to be installed or on the path, so this method is disabled.`. Download bedtools binary from official website from [here](https://bedtools.readthedocs.io/en/latest/content/installation.html#downloading-a-pre-compiled-binary) and put it in the conda bin directory. 
- GPU support for torchbiggraph. Check [this](https://github.com/facebookresearch/PyTorch-BigGraph#installation).