import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import SCMBench.traj_conserv_metrics as tm
import pathlib
import yaml



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
        "-o", "--output", dest="output", type=pathlib.Path, required=True,
        help="Path to output file (.yaml)"
    )
    return parser.parse_args()






def main(args: argparse.Namespace):
    rna=ad.read_h5ad(args.datasets[0])
    #print(args.root)
    indices = np.where(np.array(rna.obs['cell_type'] == args.root[0]))[0]
    #rna.uns["iroot"] = np.random.choice(indices)
    rna.uns["iroot"] = indices[0]

    emb=pd.read_csv(args.latents[0],header=None,index_col=0)

    emb_rna=ad.AnnData(emb.values,obs=rna.obs)

    sc.pp.pca(rna)
    sc.pp.neighbors(rna)
    
    sc.tl.diffmap(rna)
    sc.tl.dpt(rna)
    print('rna',rna)
    

    
    sc.pp.pca(emb_rna)
    sc.pp.neighbors(emb_rna)
    #emb_rna.uns["iroot"] = 0
    #sc.tl.dpt(emb_rna)
    #sc.tl.diffmap(emb_rna)
    


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