## Preprocessing

Codes for preprocessing of the datasets, including feature selection and sampling, unimodal representation, and calculation of the activity matrix. 

[Here](https://drive.google.com/drive/folders/1Ch3iPsdIIL2V_6LAk5M_TD0vhUgeHiyd?usp=sharing) we provide all the datasets in h5ad format. Before you start, please download and place the datasets under `download/`. 


We provide the jupyter notebook for feature selection and sampling, python script for unimodal representation, and calculation of the activity matrix in `preprocessing/` and the corresponding example shell script for running. You can process your own dataset follow the instructions.


`preprocessing/scATAC_source.r` was derived from Deepmaps. [Ma, A. et al. Single-cell biological network inference using a heterogeneous graph transformer. Nature Communications 14, 964 (2023).]

## Simulation data generation

We use scMultiSim for simulation data generation. More details can be found [here](https://zhanglabgt.github.io/scMultiSim/).