import os
from pathlib import Path
from scipy import io
import pandas as pd
from scanpy import AnnData

#adapated from 'https://medium.com/@daimin0514/how-to-convert-singlecellexperiment-to-anndata-8ec678b3a99e'
dir = '/home/'

def save_data_for_R(adata, save_dir, layer='counts', cell_metadata=None, gene_metadata=None):
    """
    Save single-cell data in Matrix Market format for use in R. Save counts data, cell metadata, and gene metadata in
    separate files.
    
    Args:
        - adata (AnnData): AnnData object containing single-cell data
        - save_dir (str): path to directory where output files will be saved
        - layer (str): name of the key in adata.layers that contains the counts data (default: 'counts')
        - cell_metadata (pd.DataFrame): DataFrame containing cell metadata (default: adata.obs)
        - gene_metadata (pd.DataFrame): DataFrame containing gene metadata (default: adata.var)
        
    Returns:
        None
    """
    
    # Check that input is valid
    if not isinstance(adata, AnnData):
        raise ValueError("Input 'adata' must be an AnnData object.")
    if layer not in adata.layers.keys():
        raise KeyError(f"Counts data not found in adata.layers. Please make sure counts data is stored under the key {layer}.")
    
    # Get default cell and gene metadata from adata.obs and adata.var, respectively
    if cell_metadata is None:
        cell_metadata = adata.obs.copy()
    if gene_metadata is None:
        gene_metadata = adata.var.copy()
    
    # Create directory if necessary
    Path(os.path.join(save_dir, 'data_for_R')).mkdir(parents=True, exist_ok=True)
    
    # Save counts data
    if adata.layers[layer].dtype == np.float64 or adata.layers[layer].dtype == np.float32:
        adata.layers[layer].data = adata.layers[layer].data.astype(np.int16)
    io.mmwrite(os.path.join(save_dir, 'data_for_R/counts.mtx'), adata.layers[layer])
    
    # Save cell metadata
    cell_metadata = cell_metadata.copy()
    cell_metadata['Barcode'] = cell_metadata.index
   
    if 'X_umap' in adata.obsm.keys():
        if 'UMAP1' not in cell_metadata.columns:
            cell_metadata['UMAP1'] = adata.obsm['X_umap'][:, 0]
        if 'UMAP2' not in cell_metadata.columns:
            cell_metadata['UMAP2'] = adata.obsm['X_umap'][:, 1]
    else:
        print("adata.obsm['X_umap'] does not exist.")
        
    if 'X_diffmap' in adata.obsm.keys():
        if 'DIFFMAP1' not in cell_metadata.columns:
            cell_metadata['DIFFMAP1'] = adata.obsm['X_diffmap'][:, 0]
        if 'DIFFMAP2' not in cell_metadata.columns:
            cell_metadata['DIFFMAP2'] = adata.obsm['X_diffmap'][:, 1]
    else:
        print("adata.obsm['X_diffmap'] does not exist.")

    cell_metadata.to_csv(os.path.join(save_dir, 'data_for_R/counts_cellMeta.csv'), index=None)
    
    # Save gene metadata
    gene_metadata = gene_metadata.copy()
    gene_metadata['GeneName'] = gene_metadata.index
    gene_metadata.to_csv(os.path.join(save_dir, 'data_for_R/counts_geneMeta.csv'), index=None)
    
    print(f"Single-cell data successfully saved in {save_dir}/data_for_R")


#run 
save_data_for_R(adata,
                save_dir=dir,
                layer='counts',
               )


#While loading, remember that the count matrix is to be transposed
