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
        "-l1", "--latents1", dest="latents1", type=pathlib.Path, required=True,
        nargs="+", help="Path to latent embeddings (.csv)"
    
    )
    parser.add_argument(
        "-l2", "--latents2", dest="latents2", type=pathlib.Path, required=True,
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






#datasets = [anndata.read_h5ad(item) for item in args.datasets]
#latents = [pd.read_csv(item, header=None, index_col=0).to_numpy() for item in args.latents]
def main(args: argparse.Namespace):
    rna=ad.read_h5ad(args.datasets[0])
    print('rna',rna)
    atac=ad.read_h5ad(args.datasets[1])
    print('atac',atac)
    branch_1=['HSPC','CD8 TEM_2','CD8 TEM_1','CD4 Naive','CD8 Naive','CD4 TCM','CD4 TEM','Naive B','Memory B','Intermediate B']
    # branch_1=['HSPC','CD8 TEM_2','CD8 TEM_1','CD4 Naive','CD8 Naive','CD4 TCM','CD4 TEM','Naive B','Memory B','Intermediate B','NK','CD16 Mono','CD14 Mono']
    obs_index=list(rna.obs.index)
    fake_atac_index=[x + '_atac' for x in obs_index]
    fake_atac=atac.copy()
    fake_atac.obs['index']=fake_atac_index
    fake_atac.obs.set_index('index',inplace=True)

    emb_temp=pd.read_csv(args.latents1[0],header=None,index_col=0)
    emb_rna=ad.AnnData(emb_temp.values,obs=rna.obs)
    emb_temp=pd.read_csv(args.latents2[0],header=None,index_col=0)
    emb_atac=ad.AnnData(emb_temp.values,obs=fake_atac.obs)
    emb2 = ad.concat([emb_rna, emb_atac])

    if args.root[0]=='HSPC':
        print('10x-Multiome-Pbmc10k-small branch1')
        rna=rna[rna.obs['cell_type'].isin(branch_1)]
        fake_atac=fake_atac[fake_atac.obs['cell_type'].isin(branch_1)]

    id=5
    intersection = fake_atac.var.index.intersection(rna.var.index)
    rna_intersection = rna[:,intersection].copy()
    atac_intersection = fake_atac[:,intersection].copy()
    indices = np.where(np.array(rna_intersection.obs['cell_type'] == 'HSPC'))[0]
    rna_intersection.uns["iroot"] = indices[id]
    sc.pp.pca(rna_intersection)
    sc.pp.neighbors(rna_intersection)
    sc.tl.diffmap(rna_intersection)
    sc.tl.dpt(rna_intersection)
    # print(rna_intersection.uns['diffmap_evals'].shape)
    indices = np.where(np.array(atac_intersection.obs['cell_type'] == 'HSPC'))[0]
    atac_intersection.uns["iroot"] = indices[id]
    sc.pp.pca(atac_intersection)
    sc.pp.neighbors(atac_intersection)
    sc.tl.diffmap(atac_intersection)
    sc.tl.dpt(atac_intersection)
    print(atac_intersection.uns['diffmap_evals'].shape)
    

    rna2 = ad.concat([rna_intersection, atac_intersection])
    rna2.uns["iroot"] = indices[id]

    if args.root[0]=='HSPC':
        print('10x-Multiome-Pbmc10k-small branch1')
        emb2=emb2[emb2.obs['cell_type'].isin(branch_1)]
    sc.pp.pca(emb2)
    sc.pp.neighbors(emb2)

    rna2.write_h5ad(str(args.datasets[0])[:-5]+'_fake.h5ad')
    emb2.write_h5ad(str(args.latents1[0])[:-4]+'_fake.h5ad')
    #emb_rna.uns["iroot"] = 0
    #sc.tl.dpt(emb_rna)
    #sc.tl.diffmap(emb_rna)
    

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
    
    
    #print(score)


if __name__ == "__main__":
    main(parse_args())