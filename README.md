# SCMBench: Benchmarking Domain-specific and Foundation Models for Single-cell Multi-omics Data Integration



## Directory structure

```
.
├── SCMBench                # Main Python package
├── data                    # Data files
├── evaluation              # Method evaluation pipelines
├── methods                 # Tools included in the benchmarking
├── environments            # Reproducible Python or R environment
├── pyproject.toml          # Python package metadata
├── LICENSE
└── README.md
```

## Installation

### 1. Clone our repository

```
git clone https://github.com/Susanxuan/SCMBench.git
```

### 2. Set up for each method

- We have summarized Python/R version and packages needed for each method in [environments](environments), users can `bash [env]_installation.sh` files to install requirements for each method accordingly.

| Virtual Envs         | Python/R Version | Methods                                                                             |
| -------------------- | ---------------- | ----------------------------------------------------------------------------------- |
| env1_installation.sh | python=3.8       | cobolt;GLUE;Harmony;MMD-MA;MOFA;Pamona;PCA;<br />scJoint;scMDC;scMoMaT;UCE;UnionCom |
| env2_installation.sh | python=3.11      | Geneformer;scFoundation;scVI;TotalVI                                                |
| env3_installation.sh | python=3.10      | scGPT                                                                               |
| env4_installation.sh | R=4.3.3          | bindSC;Deepmaps;iNMF;liger;Seurat4;Seurat5   


## Quick start
1. Data Preprocessing
Follow [this](/data/README.md).

2. Run python scrips. E.g., scJoint.\:
```bash
cd methods/scJoint
mkdir ../results/scJoint-output
python run_scJoint.py --input-rna ../data/download/10x-Multiome-Pbmc10k-small-RNA.h5ad --input-atac  10x-Multiome-Pbmc10k-small-ATAC.h5ad  --output-rna ../results/scJoint-output/rna.csv --output-atac ../results/scJoint-output/atac.csv  -r ../results/scJoint-output/run_info.yaml
```

3. Run R script. E.g., Deepmaps.\
Tested approach: install R and related packages with in conda (e.g. [this](https://stackoverflow.com/questions/70410968/is-it-possible-to-install-r-in-miniconda)).
```bash
cd methods/Deepmaps
mkdir ../results/Deepmaps-output
Rscript run_Deepmaps.R --input-rna ../data/download/10x-Multiome-Pbmc10k-small-RNA.h5ad --input-atac  10x-Multiome-Pbmc10k-small-ATAC.h5ad  --output-rna ../results/Deepmaps-output/rna.csv --output-atac ../results/Deepmaps-output/atac.csv  -r ../results/Deepmaps-output/run_info.yaml
```

## Data preprocessing

Preprocess all the data to the same `.h5ad` format with the same keys following [here](data/README.md).

For new datasets we include, upload preprocessing scripts and update the preprocessed data/link [here](data/README.md) (store the preprocessed data somewhere we can directly download).

Current included datasets:
- ├─ 10x-Multiome-Pbmc10k\
  ├─ Chen-2019\
  ├─ Ma-2020\
  ├─ Muto-2021\
  ├─ Yao-2021\
  └─ Triple


## Algorithms

Note: Most tools can be directly used by installing it individually via `pip install` if they provide that option. But for scJoint, and UCE, we provide edited package in [methods/](methods/) used for their corresponding scripts `run_[METHOD].py`.

Current included algorithms:

### Statistical-Based:

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
- [Harmony](https://www.nature.com/articles/s41592-019-0619-0)

### Deep Learning-based:

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

### Foundation Models:
- [Geneformer](https://www.nature.com/articles/s41586-023-06139-9)
- [scFoundation](https://www.nature.com/articles/s41592-024-02305-7)
- [scGPT](https://www.nature.com/articles/s41592-024-02201-0)
- [UCE](https://www.biorxiv.org/content/10.1101/2023.11.28.568918v1)

## Downstream tasks

For downstream evaluation, please read [evaluation/README.md](evaluation/README.md).

Currently included downstream tasks:
- Multi-Omic Integration Accuracy (MAP, MNI, ASW, ARI)
- Bio-Conservation: 
  - Biomarker Detection (JSI)
  - DARs Detection (JSI)
  - Enriched Motifs Detection (JSI)
  - Trajectory Conservation
- Batch effect & detecting over-correction (iLISI, kBET, Graph connectivity, Batch-ASW)

## Solution to potential install issues

- Bedtools. `NotImplementedError: "intersectBed" does not appear to be installed or on the path, so this method is disabled.`. Download bedtools binary from official website from [here](https://bedtools.readthedocs.io/en/latest/content/installation.html#downloading-a-pre-compiled-binary) and put it in the conda bin directory. 
- GPU support for torchbiggraph. Check [this](https://github.com/facebookresearch/PyTorch-BigGraph#installation).
