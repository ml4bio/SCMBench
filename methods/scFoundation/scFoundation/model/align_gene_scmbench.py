import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
from typing import List, Union, Optional
from scipy import sparse
import argparse

def align_expression_data(
    expr_data: Union[np.ndarray, sparse.spmatrix],
    gene_names: List[str],
    reference_genes: List[str]
) -> np.ndarray:
    """
    Align expression data with reference gene list
    
    Parameters:
    -----------
    expr_data : Union[np.ndarray, sparse.spmatrix]
        Expression matrix to align
    gene_names : List[str]
        List of gene names corresponding to expr_data columns
    reference_genes : List[str]
        Ordered list of reference genes to align to
        
    Returns:
    --------
    numpy.ndarray
        Aligned expression matrix
    """
    # Convert sparse matrix to dense if needed
    if sparse.issparse(expr_data):
        expr_data = expr_data.toarray()
        
    # Create DataFrame from expression data
    expr_df = pd.DataFrame(expr_data, columns=gene_names)
    
    # Create new DataFrame with reference genes
    new_expr_df = pd.DataFrame(0, 
                             index=expr_df.index,
                             columns=reference_genes)
    
    # Find common genes and fill expression values while preserving gene order
    common_genes = list(set(reference_genes) & set(expr_df.columns))
    if not common_genes:
        raise ValueError("No common genes found between expression data and reference genes")
    
    new_expr_df[common_genes] = expr_df[common_genes]
            
    return new_expr_df.values

def process_and_align_data(args):
    """
    Process data splits and align genes with reference
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments containing:
        - reference_path: path to reference gene list
        - data_path: path to input data
        - data_name: dataset name
        - file_suffix: suffix of the h5ad file ('ACTIVE.h5ad' or 'RNA.h5ad')
        - use_hvg: whether to use highly variable genes
        - n_top_features: number of top variable features to select (if use_hvg is True)
        - output_dir: directory for output files
        
    Returns:
    --------
    numpy.ndarray
        Aligned expression data
    """
    # Input validation
    if not os.path.exists(args.reference_path):
        raise FileNotFoundError(f"Reference file not found: {args.reference_path}")
    
    load_path = os.path.join(args.data_path, args.data_name, args.data_name + '-' + args.file_suffix)
    if not os.path.exists(load_path):
        raise FileNotFoundError(f"Data file not found: {load_path}")
        
    # Load the reference gene list
    print("Reading reference gene list...")
    try:
        ref_df = pd.read_csv(args.reference_path, sep='\t')
        reference_genes = ref_df['gene_name'].tolist()
    except Exception as e:
        raise ValueError(f"Error reading reference gene list: {str(e)}")
    
    # Load and process data
    print(f'Loading {load_path}')
    try:
        adata = ad.read_h5ad(load_path)
    except Exception as e:
        raise ValueError(f"Error reading expression data: {str(e)}")
    
    # Log statistics of the raw data
    print(f'Raw data shape: {adata.shape}')
    
    # Count non-zero genes in raw data
    if sparse.issparse(adata.X):
        non_zero_raw = (adata.X != 0).sum()
        if hasattr(non_zero_raw, 'A'):  # For CSR matrix
            non_zero_raw = non_zero_raw.A[0]
        print(f'Number of non-zero values in raw data: {non_zero_raw}')
    else:
        non_zero_raw = np.count_nonzero(adata.X)
        print(f'Number of non-zero values in raw data: {non_zero_raw}')
    
    # Feature selection if use_hvg is True
    if args.use_hvg:
        print(f'Selecting top {args.n_top_features} highly variable genes')
        
        # Make a copy and ensure we're working with dense array for preprocessing
        org_adata = adata.copy()
        if sparse.issparse(org_adata.X):
            org_adata.X = org_adata.X.toarray()
        
        # Log transform and select highly variable genes
        sc.pp.log1p(org_adata)
        sc.pp.highly_variable_genes(
            org_adata,
            n_top_genes=args.n_top_features,
            flavor='seurat',
            span=0.3
        )
        
        # Apply feature selection mask
        hvg_mask = org_adata.var['highly_variable']
        hvg_genes = org_adata.var_names[hvg_mask].tolist()
        
        print(f'Selected {len(hvg_genes)} highly variable genes')
        
        # Filter the original data to keep only HVG
        adata = adata[:, hvg_genes]
        
        # Add suffix to output filename
        output_suffix = f'_hvg{args.n_top_features}'
    else:
        print('Skipping highly variable gene selection')
        output_suffix = ''
    
    print('\nProcessing and aligning splits...')
    try:
        aligned_data = align_expression_data(
            adata.X,
            adata.var_names.tolist(),
            reference_genes
        )
    except Exception as e:
        raise ValueError(f"Error in alignment: {str(e)}")
    
    # Create and save DataFrame
    aligned_df = pd.DataFrame(aligned_data, columns=reference_genes)
    aligned_df[aligned_df.abs() < 1e-6] = 0
    aligned_df = aligned_df.round(6)
    
    # Save results
    output_path = os.path.join(
        args.output_dir, 
        f'{args.data_name}{output_suffix}_aligned.csv'
    )
    
    try:
        aligned_df.to_csv(output_path, index=False)
        
        # Also save statistics to a log file
        log_path = os.path.join(
            args.output_dir, 
            f'{args.data_name}{output_suffix}_stats.log'
        )
        with open(log_path, 'w') as log_file:
            log_file.write(f"Dataset: {args.data_name}\n")
            log_file.write(f"Use HVG: {args.use_hvg}\n")
            if args.use_hvg:
                log_file.write(f"Number of top HVG: {args.n_top_features}\n")
            log_file.write(f"Raw data shape: {adata.shape}\n")
            log_file.write(f"Aligned data shape: {aligned_df.shape}\n")
            log_file.write(f"Number of non-zero genes: {(aligned_df != 0).any().sum()}\n")
            log_file.write(f"Average non-zero genes per cell: {(aligned_df != 0).sum(axis=1).mean():.2f}\n")
            log_file.write(f"Min non-zero genes per cell: {(aligned_df != 0).sum(axis=1).min()}\n")
            log_file.write(f"Max non-zero genes per cell: {(aligned_df != 0).sum(axis=1).max()}\n")
            log_file.write(f"Percentage of values that are non-zero: {(aligned_df != 0).sum().sum() / (aligned_df.shape[0] * aligned_df.shape[1]) * 100:.2f}%\n")
        print(f"Statistics saved to {log_path}")
        
    except Exception as e:
        raise IOError(f"Error saving output files: {str(e)}")
    
    # Print statistics
    print(f'Original shape: {aligned_df.shape}')
    non_zero_genes = (aligned_df != 0).any().sum()
    print(f'Number of non-zero genes: {non_zero_genes}')
    
    # Calculate non-zero genes per cell
    non_zero_per_cell = (aligned_df != 0).sum(axis=1)
    print(f'Average non-zero genes per cell: {non_zero_per_cell.mean():.2f}')
    print(f'Min non-zero genes per cell: {non_zero_per_cell.min()}')
    print(f'Max non-zero genes per cell: {non_zero_per_cell.max()}')
    
    return aligned_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', default='/mnt/nas/user/yixuan/SCMBench/data/download/', type=str,
                       help='Location where data is stored')
    parser.add_argument('--data_name', type=str, required=True,
                       help='Dataset name to load (e.g., Muto-2021-small)')
    parser.add_argument('--file_suffix', type=str, default='ACTIVE.h5ad',
                       help='File suffix (e.g., ACTIVE.h5ad or RNA.h5ad)')
    parser.add_argument('--reference_path', type=str, required=True,
                       help='Path to reference gene list TSV file')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Directory to save aligned data')
    parser.add_argument('--use_hvg', action='store_true',
                       help='Whether to use highly variable genes')
    parser.add_argument('--n_top_features', type=int, default=1000,
                       help='Number of top variable features to select (only used if use_hvg is set)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process and align data
    aligned_data = process_and_align_data(args)
    
    print('\nProcessing complete!')

if __name__ == '__main__':
    main()

'''
python align_gene_scmbench.py \
  --data_name Muto-2021-batch-1-small\
  --file_suffix ACTIVE.h5ad \
  --reference_path OS_scRNA_gene_index.19264.tsv \
  --output_dir /mnt/nas/user/yixuan/scFoundation/model/scmdata/aligned_data \
  --use_hvg --n_top_features 5000

CUDA_VISIBLE_DEVICES=1 python get_embedding.py --task_name scmbench --input_type singlecell --output_type cell --pool_type all --tgthighres a5 \
    --data_path /mnt/nas/user/yixuan/scFoundation/model/scmdata/aligned_data/Muto-2021-batch-1-small_aligned.csv \
    --save_path /mnt/nas/user/yixuan/scFoundation/model/scmdata/emb \
    --pre_normalized T --version rde

'''
