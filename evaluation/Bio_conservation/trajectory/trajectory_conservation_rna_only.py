import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import pathlib
import yaml
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
scmbench_dir = os.path.abspath(os.path.join(parent_dir,'..'))
print('scmbench_dir',scmbench_dir)
sys.path.append(scmbench_dir)

import SCMBench.traj_conserv_metrics as tm


def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Compute trajectory conservation metrics for paired samples"
    )
    parser.add_argument(
        "-d", "--datasets", dest="datasets", type=pathlib.Path, required=True,
        nargs="+", help="Path to datasets (.h5ad)"
    )
    parser.add_argument(
        "-l", "--latents", dest="latents", type=pathlib.Path, required=True,
        nargs="+", help="Path to latent embeddings (.csv)"
    
    )
    parser.add_argument(
        "-r", "--root", dest="root", type=str, required=True,
        nargs="+", help="root cell type"
    
    )
    parser.add_argument(
        "-b", "--branch", dest="branch", action="store_true", 
        help="whether specify branch"
    
    )
    parser.add_argument(
        "-bn", "--branch_name", dest="branch_name", type=str,
        nargs="+", help="cell types name on the branch"
    
    )
    parser.add_argument(
        "-c", "--comb_data", dest="comb_data", type=pathlib.Path,
        help="whether the combined data file (with calculated traj info) is provided, if provided skip calculation from input datasets"
    
    )
    parser.add_argument(
        "-o", "--output", dest="output", type=pathlib.Path, required=True,
        help="Path to output file (.yaml)"
    )
    return parser.parse_args()






def main(args: argparse.Namespace):
    print("[1/3] Reading in data...")
    rna=ad.read_h5ad(args.datasets[0])
    
    # define the root indices by your self
    indices = np.where(np.array(rna.obs['cell_type'] == args.root[0]))[0]
    #rna.uns["iroot"] = np.random.choice(indices)
    rna.uns["iroot"] = indices[0] #scjoiny idx 1

    emb=pd.read_csv(args.latents[0],header=None,index_col=0)
    emb_rna=ad.AnnData(emb.values,obs=rna.obs)

    if args.comb_data:
        if str(args.comb_data).endswith('.h5ad'):
            print('Combined data file (with calculated traj info in h5ad) is provided, skip redundant calculation')
            rna=ad.read_h5ad(args.comb_data)
        else: 
            print('Combined data file should be in h5ad format')
    else:
        print('Combined data file (with calculated traj info in h5ad) is not provided, calculating traj info GT')

        if args.branch:
            print('Specify branch:',args.branch_name)
            rna=rna[rna.obs['cell_type'].isin(args.branch_name)]

        sc.pp.pca(rna)
        sc.pp.neighbors(rna)
        
        sc.tl.diffmap(rna)
        sc.tl.dpt(rna)
        rna_path = os.path.join(os.path.dirname(args.output),'raw_combine_traj_rna_only_branch1_idx0.h5ad')
        rna.write_h5ad(rna_path)

    if args.branch:
        print('10x-Multiome-Pbmc10k-small branch1')
        emb_rna=emb_rna[emb_rna.obs['cell_type'].isin(args.branch_name)]
    
    sc.pp.pca(emb_rna)
    sc.pp.neighbors(emb_rna)
    
    # save emb for visualization
    # emb_rna.write_h5ad(str(args.latents[0])[:-4]+'_traj_rna_emb.h5ad')
    


    score = tm.trajectory_conservation(
            adata_pre=rna,
            adata_post=emb_rna,
            label_key="cell_type",
            pseudotime_key="dpt_pseudotime",
        )

    metrics = {
        "trajectory_conservation_score": score,
    }
    
    # round results to .4f
    for k, v in metrics.items():
        metrics[k] = round(float(v), 4)
    print("[3/3] Saving results...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        yaml.dump(metrics, f)
    
    

if __name__ == "__main__":
    main(parse_args())