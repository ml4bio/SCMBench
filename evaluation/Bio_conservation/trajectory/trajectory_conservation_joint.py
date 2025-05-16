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
        help="root cell type"
    
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
    atac=ad.read_h5ad(args.datasets[1])
    obs_index=list(rna.obs.index)
    activity_index=[x + '_atac' for x in obs_index]
    activity=atac.copy()
    activity.obs['index']=activity_index
    activity.obs.set_index('index',inplace=True)

    emb_temp=pd.read_csv(args.latents[0],header=None,index_col=0)
    emb_rna=ad.AnnData(emb_temp.values,obs=rna.obs)
    emb_temp=pd.read_csv(args.latents[1],header=None,index_col=0)
    emb_atac=ad.AnnData(emb_temp.values,obs=activity.obs)
    emb2 = ad.concat([emb_rna, emb_atac])
    
    if args.comb_data:
        if str(args.comb_data).endswith('.h5ad'):
            print('Combined data file (with calculated traj info in h5ad) is provided, skip redundant calculation')
            rna2=ad.read_h5ad(args.comb_data)
        else: 
            print('Combined data file should be in h5ad format')
    else:
        print('Combined data file (with calculated traj info in h5ad) is not provided, calculating traj info GT')

        if args.branch:
            print('Specify branch:',args.branch_name)
            rna=rna[rna.obs['cell_type'].isin(args.branch_name)]
            activity=activity[activity.obs['cell_type'].isin(args.branch_name)]

        id=5
        intersection = activity.var.index.intersection(rna.var.index)
        rna_intersection = rna[:,intersection].copy()
        atac_intersection = activity[:,intersection].copy()
        indices = np.where(np.array(rna_intersection.obs['cell_type'] == 'HSPC'))[0]
        rna_intersection.uns["iroot"] = indices[id]
        sc.pp.pca(rna_intersection)
        sc.pp.neighbors(rna_intersection)
        sc.tl.diffmap(rna_intersection)
        sc.tl.dpt(rna_intersection)
        indices = np.where(np.array(atac_intersection.obs['cell_type'] == 'HSPC'))[0]
        atac_intersection.uns["iroot"] = indices[id]
        sc.pp.pca(atac_intersection)
        sc.pp.neighbors(atac_intersection)
        sc.tl.diffmap(atac_intersection)
        sc.tl.dpt(atac_intersection)
        print(atac_intersection.uns['diffmap_evals'].shape)
        

        rna2 = ad.concat([rna_intersection, atac_intersection])
        rna2.uns["iroot"] = indices[id]
        rna2_path = os.path.join(os.path.dirname(args.output),'raw_combine_traj.h5ad')
        rna2.write_h5ad(rna2_path)
        print('Combined data file (with calculated traj info in h5ad) saved to ',rna2_path)
        print('idx used',id)

    print("[2/3] Calculating...")

    if args.branch:
        print('10x-Multiome-Pbmc10k-small branch1')
        emb2=emb2[emb2.obs['cell_type'].isin(args.branch_name)]
    sc.pp.pca(emb2)
    sc.pp.neighbors(emb2)

    # save emb for visualization
    # emb2.write_h5ad(str(args.latents[0])[:-7]+'combo_emb.h5ad')
    

    score = tm.trajectory_conservation(
            adata_pre=rna2,
            adata_post=emb2,
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