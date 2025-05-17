"""
Run scMDC
"""

import argparse
from time import time
import math, os
from sklearn import metrics
from sklearn.cluster import KMeans
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.nn import Parameter
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

sys.path.append('scMDC/src')
from scMDC import scMultiCluster
import numpy as np
import pandas as pd
import collections
import h5py
import scanpy as sc
from preprocess import read_dataset, normalize, clr_normalize_each_cell
from sklearn.preprocessing import LabelEncoder
from utils import *
import yaml
import pathlib


def parse_args() -> argparse.Namespace:
    r"""
    Parse command line arguments
    """
    # setting the hyper parameters
    
    parser = argparse.ArgumentParser(description='train',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_clusters', default=27, type=int)
    parser.add_argument('--cutoff', default=0.5, type=float, help='Start to train combined layer after what ratio of epoch')
    parser.add_argument('--batch_size', default=256, type=int)
    parser.add_argument('--input-rna', default='/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/Chen-2019/Chen-2019-RNA.h5ad')
    parser.add_argument('--input-atac', default='/mnt/nas/user/yixuan/Multiomics-benchmark-main/data/download/Chen-2019/Chen-2019-ATAC.h5ad')
    parser.add_argument('--maxiter', default=5000, type=int)
    parser.add_argument('--pretrain_epochs', default=400, type=int)
    parser.add_argument('--gamma', default=.1, type=float,
                        help='coefficient of clustering loss')
    parser.add_argument('--tau', default=1., type=float,
                        help='fuzziness of clustering loss')                    
    parser.add_argument('--phi1', default=0.001, type=float,
                        help='coefficient of KL loss in pretraining stage')
    parser.add_argument('--phi2', default=0.001, type=float,
                        help='coefficient of KL loss in clustering stage')
    parser.add_argument('--update_interval', default=1, type=int)
    parser.add_argument('--tol', default=0.001, type=float)
    parser.add_argument('--lr', default=1., type=float)
    parser.add_argument('--ae_weights', default=None)
    parser.add_argument('--save_dir', default='results/')
    parser.add_argument('--ae_weight_file', default='AE_weights_1.pth.tar')
    parser.add_argument('--resolution', default=0.2, type=float)
    parser.add_argument('--n_neighbors', default=30, type=int)
    parser.add_argument('--embedding_file', action='store_true', default=False)
    parser.add_argument('--prediction_file', action='store_true', default=False)
    parser.add_argument('-el','--encodeLayer', nargs='+', default=[256,64,32,16])
    parser.add_argument('-dl1','--decodeLayer1', nargs='+', default=[16,64,256])
    parser.add_argument('-dl2','--decodeLayer2', nargs='+', default=[16,20])
    parser.add_argument('--sigma1', default=2.5, type=float)
    parser.add_argument('--sigma2', default=1.5, type=float)
    parser.add_argument('--f1', default=1000, type=float, help='Number of mRNA after feature selection')
    parser.add_argument('--f2', default=2000, type=float, help='Number of ADT/ATAC after feature selection')
    parser.add_argument('--filter1', action='store_true', default=False, help='Do mRNA selection')
    parser.add_argument('--filter2', action='store_true', default=False, help='Do ADT/ATAC selection')
    parser.add_argument('--run', default=1, type=int)
    parser.add_argument('--device', default='cuda')
    parser.add_argument(
        "--output-rna", dest="output_rna", type=pathlib.Path, required=True,
        help="Path of output RNA latent file (.csv)"
    )
    parser.add_argument(
        "--output-atac", dest="output_atac", type=pathlib.Path, required=True,
        help="Path of output ATAC latent file (.csv)"
    )
    parser.add_argument(
        "-r", "--run-info", dest="run_info", type=pathlib.Path, required=True,
        help="Path of output run info file (.yaml)"
    )
    return parser.parse_args()

def main(args: argparse.Namespace) -> None:
    r"""
    Main function
    """

    print("[1/4] Reading data...")
    X1_adata = sc.read_h5ad(args.data_file_1)
    X2_adata = sc.read_h5ad(args.data_file_2)
    x1 = X1_adata.X.toarray()
    x2 = X2_adata.X.toarray()
    y_str = X1_adata.obs['cell_type']
    y = LabelEncoder().fit_transform(y_str)
    
    print("[2/4] Preprocessing...")
    #Gene filter
    if args.filter1:
        importantGenes = geneSelection(x1, n=args.f1, plot=False)
        x1 = x1[:, importantGenes]
    if args.filter2:
        importantGenes = geneSelection(x2, n=args.f2, plot=False)
        x2 = x2[:, importantGenes]
        
    # preprocessing scRNA-seq read counts matrix
    adata1 = sc.AnnData(x1)
    adata1.obs['Group'] = y

    adata1 = read_dataset(adata1,
                     transpose=False,
                     test_split=False,
                     copy=True)

    adata1 = normalize(adata1,
                      size_factors=True,
                      normalize_input=True,
                      logtrans_input=True)
    
    adata2 = sc.AnnData(x2)
    adata2.obs['Group'] = y
    adata2 = read_dataset(adata2,
                     transpose=False,
                     test_split=False,
                     copy=True)
    
    adata2 = normalize(adata2,
                      size_factors=True,
                      normalize_input=True,
                      logtrans_input=True)

    #adata2 = clr_normalize_each_cell(adata2)

    input_size1 = adata1.n_vars
    input_size2 = adata2.n_vars
    
    # print(args)
    print("[3/4] Training scMDC...")
    encodeLayer = list(map(int, args.encodeLayer))
    decodeLayer1 = list(map(int, args.decodeLayer1))
    decodeLayer2 = list(map(int, args.decodeLayer2))
    
    model = scMultiCluster(input_dim1=input_size1, input_dim2=input_size2, tau=args.tau,
                        encodeLayer=encodeLayer, decodeLayer1=decodeLayer1, decodeLayer2=decodeLayer2,
                        activation='elu', sigma1=args.sigma1, sigma2=args.sigma2, gamma=args.gamma, 
                        cutoff = args.cutoff, phi1=args.phi1, phi2=args.phi2, device=args.device).to(args.device)
    
    # print(str(model))
    
    if not os.path.exists(args.save_dir):
            os.makedirs(args.save_dir)
            
    start_time = time()
    if args.ae_weights is None:
        model.pretrain_autoencoder(X1=adata1.X, X_raw1=adata1.raw.X, sf1=adata1.obs.size_factors, 
                X2=adata2.X, X_raw2=adata2.raw.X, sf2=adata2.obs.size_factors, batch_size=args.batch_size, 
                epochs=args.pretrain_epochs, ae_weights=args.ae_weight_file)
    else:
        if os.path.isfile(args.ae_weights):
            print("==> loading checkpoint '{}'".format(args.ae_weights))
            checkpoint = torch.load(args.ae_weights)
            model.load_state_dict(checkpoint['ae_state_dict'])
        else:
            print("==> no checkpoint found at '{}'".format(args.ae_weights))
            raise ValueError
    
    # print('Pretraining time: %d seconds.' % int(time() - t0))
    
    #get k
    latent = model.encodeBatch(torch.tensor(adata1.X).to(args.device), torch.tensor(adata2.X).to(args.device))
    latent = latent.cpu().numpy()
    if args.n_clusters == -1:
       n_clusters = GetCluster(latent, res=args.resolution, n=args.n_neighbors)
    else:
       print("n_cluster is defined as " + str(args.n_clusters))
       n_clusters = args.n_clusters

    y_pred, _ = model.fit(X1=adata1.X, X_raw1=adata1.raw.X, sf1=adata1.obs.size_factors, 
        X2=adata2.X, X_raw2=adata2.raw.X, sf2=adata2.obs.size_factors, y=y,
        n_clusters=n_clusters, batch_size=args.batch_size, num_epochs=args.maxiter, 
        update_interval=args.update_interval, tol=args.tol, lr=args.lr, save_dir=args.save_dir)
    elapsed_time = time() - start_time
    # print('Total time: %d seconds.' % int(time() - t0))

    print("[4/4] Saving results...")
    if args.prediction_file:
       y_pred_ = best_map(y, y_pred) - 1
    #    np.savetxt(args.save_dir + "/" + str(args.run) + "_pred.csv", y_pred_, delimiter=",")
    
    if args.embedding_file:
       final_latent = model.encodeBatch(torch.tensor(adata1.X).to(args.device), torch.tensor(adata2.X).to(args.device))
       final_latent = final_latent.cpu().numpy()

    rna_latent = final_latent
    atac_latent = final_latent
    #    np.savetxt(args.save_dir + "/" + str(args.run) + "_embedding.csv", final_latent, delimiter=",")

    args.output_rna.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rna_latent, index=X1_adata.obs_names).to_csv(args.output_rna, header=False)
    args.output_atac.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(atac_latent, index=X2_adata.obs_names).to_csv(args.output_atac, header=False)
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    with args.run_info.open("w") as f:
        yaml.dump({
            "cmd": " ".join(sys.argv),
            "args": vars(args),
            "time": elapsed_time,
            "n_cells": x1.shape[0] + x2.shape[0]
        }, f)


if __name__ == "__main__":
    main(parse_args())