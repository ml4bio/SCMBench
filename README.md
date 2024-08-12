# Dev Code for Multiomics Benchmark

Code adapted from [GLUE](https://github.com/gao-lab/GLUE). Major refactoring needed!

Some usage details can also be found in their [docs](https://scglue.readthedocs.io)



## Directory structure

```
.
├── SCMBench                  # Main Python package
├── data                    # Data files
├── evaluation              # Method evaluation pipelines
├── experiments             # Experiments and case studies
├── tests                   # Unit tests for the Python package
├── docs                    # Documentation files
├── custom                  # Customized third-party packages
├── packrat                 # Reproducible R environment via packrat
├── env.yaml                # Reproducible Python environment via conda
├── pyproject.toml          # Python package metadata
├── LICENSE
└── README.md
```

Can ignore the `snakemake` for now.

## Quick start
1. Data Preprocessing
Follow [this](https://SCMBench.readthedocs.io/en/latest/preprocessing.html).

2. Run GLUE:
```bash
python setup.py develop
cd evaluation/workflow/scripts
python run_GLUE.py --input-rna /home/co-zong1/rds/hpc-work/multiomics/Multiomics-benchmark/data/download/Chen-2019/Chen-2019-RNA.h5ad --input-atac  /home/co-zong1/rds/hpc-work/multiomics/Multiomics-benchmark/data/download/Chen-2019/Chen-2019-ATAC-preprocessed.h5ad -p /home/co-zong1/rds/hpc-work/multiomics/Multiomics-benchmark/data/download/Chen-2019/guidance.graphml.gz --train-dir ./glue-output --output-rna ./glue-output/rna.csv --output-atac ./glue-output/atac.csv --output-feature ./glue-output/features.csv -r glue-output/run_info.yaml
```

3. Run R script. E.g., Harmony.\
Tested approach: install R and related packages with in conda (e.g. [this](https://stackoverflow.com/questions/70410968/is-it-possible-to-install-r-in-miniconda)).
```bash
cd evaluation/workflow/scripts
mkdir harmony-output
Rscript run_Harmony.R --input-rna /home/co-zong1/rds/hpc-work/multiomics/Multiomics-benchmark/data/download/Chen-2019/Chen-2019-RNA.h5ad --input-atac  /home/co-zong1/rds/hpc-work/multiomics/Multiomics-benchmark/data/download/Chen-2019/Chen-2019-ATAC-preprocessed.h5ad --output-rna ./harmony-output/rna.csv --output-atac ./harmony-output/atac.csv --run-info harmony-output/run_info.yaml
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

For newly added downstream tasks, follow the code fashion of `run_[METHOD].py` in [evaluation/workflow/scripts/](evaluation/workflow/scripts).

Note: you can directly import new algorithms via `pip install` if they provide that option. But do remember to record the version and update `env.yaml`.

Current included algorithms:
- [GLUE](https://www.nature.com/articles/s41587-022-01284-4)
- [CLUE](https://openreview.net/forum?id=Tfb73TeKnJ-)
- [bindSC](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02679-x)
- [Harmony](https://www.nature.com/articles/s41592-019-0619-0)
- [Seurat v3](https://www.cell.com/cell/pdf/S0092-8674(19)30559-8.pdf)
- [iNMF](https://www.nature.com/articles/s41587-021-00867-x)
- [LIGER](https://www.cell.com/cell/pdf/S0092-8674(19)30504-5.pdf)
- [MMD_MA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8496402/)
- [Pamona](https://academic.oup.com/bioinformatics/article/38/1/211/6353029)
- [UnionCom](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i48/5870490)
- [SIMBA](https://www.nature.com/articles/s41592-023-01899-8)

## Downstream tasks

For newly added downstream tasks, follow the code fashion in [Experiments/](experiments).

Currently included downstream tasks:
- Cell type
- Batch effect & detecting over-correction
- Regulatory inference

## Solution to potential install issues

- Bedtools. `NotImplementedError: "intersectBed" does not appear to be installed or on the path, so this method is disabled.`. Download bedtools binary from official website from [here](https://bedtools.readthedocs.io/en/latest/content/installation.html#downloading-a-pre-compiled-binary) and put it in the conda bin directory. 
- GPU support for torchbiggraph. Check [this](https://github.com/facebookresearch/PyTorch-BigGraph#installation).