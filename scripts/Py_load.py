# Script to load files saved from seurat into python as an Anndata object.

import scanpy as sc
import pandas as pd
from scipy import io, sparse

save_dir = ''


def loadDataSavedFromSeurat(path: str, prefix: str = "RNA"):
    """
    Load data from a specified path.

    Args:
        path (str): Path to the directory containing the input files.
        prefix (str, optional): Prefix for the input file names (default: 'RNA').

    Returns:
        Anndata object with loaded counts, barcodes, genes, cell meta data and umap coordinates.
    """
    
    # Load the count matrix from the provided path
    try:
        counts = io.mmread(f"{path}/{prefix}_counts.mtx")
        counts = sparse.csr_matrix(counts)
    except FileNotFoundError:
        raise FileNotFoundError("Count matrix file not found. Please check the path and prefix.")
    
    # Load the barcodes and genes information
    try:
        barcodes = pd.read_csv(f"{path}/{prefix}_barcodes.csv")
        genes = pd.read_csv(f"{path}/{prefix}_genes.csv")
    except FileNotFoundError:
        raise FileNotFoundError("Barcodes or genes file not found. Please check the path and prefix.")
    
    # Create an AnnData object using the count matrix
    adata = sc.AnnData(counts.T)
    
    
    # Set the observation names to the barcodes
    adata.obs_names = barcodes['Barcode'].values
    
    # Set the variable names to the gene names
    adata.var_names = genes['Gene'].values
    
    # Load the cell metadata and set it to the observation annotations in adata object
    try:
        cellMeta = pd.read_csv(f"{path}/{prefix}_cellMeta.csv", index_col=0)
        adata.obs = cellMeta
    except FileNotFoundError:
        raise FileNotFoundError("Cell metadata file not found. Please check the path and prefix.")
        
    # Get UMAP coordinates from cell metadata and set them to the observation embeddings in adata object
    if 'UMAP_1' in cellMeta.columns and 'UMAP_2' in cellMeta.columns:
        adata.obsm['X_umap'] = adata.obs.loc[:, ['UMAP_1', 'UMAP_2']].values.copy()
    else:
        print("UMAP coordinates not found. Skipping UMAP coordinates...")

    # Delete the variables that are no longer needed to free up memory
    del counts, barcodes, genes, cellMeta

    return adata

adata=loadDataSavedFromSeurat(path=save_dir+'/data4py',
                        prefix='RNA'
                       )
