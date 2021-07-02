import pandas as pd
import scipy
import numpy as np
import scanpy as sc
import scanpy.external as sce
from matplotlib import pyplot as plt

def batch_process(adata_batch, remove_extreme_cells=True, scaling=True):
    if scipy.sparse.issparse(adata_batch.X):
        adata_batch.X = adata_batch.X.todense()
                
        adata_batch.X = np.asarray(adata_batch.X)
    
    adata_batch.raw = adata_batch
    
    if 'batch' in adata_batch.obs:
        adata_batch.obs['donor'] = [i[0] for i in adata_batch.obs['batch']]
        adata_batch.obs['day'] = [str(int(i[-1]) - 1) for i in adata_batch.obs['batch']]
        batches = list(dict.fromkeys(adata_batch.obs['batch'].values))
        
        if remove_extreme_cells:
            # Remove top and bottom 0.3% of cells
            keep_cells = []

            for batch in batches:
                idx = adata_batch[adata_batch.obs['batch'] == batch].obs_names 
                counts = adata_batch[adata_batch.obs['batch'] == batch].X.sum(1)
                keep_cells += list(idx[(counts < np.percentile(counts, 99.5)) & (counts > np.percentile(counts, 0.3))])

            print('Removing {} cells for extreme expression values'.format(len(adata_batch) - len(keep_cells)))

            adata_batch = adata_batch[keep_cells]
        
        if scaling:
            # Apply scaling
            dict_counts_before = {batch: adata_batch[adata_batch.obs['batch'] == batch].X.sum(1) for batch in batches}

            n_reads_per_cell = {batch: adata_batch[adata_batch.obs['batch'] == batch].X.sum() / 
                                len(adata_batch[adata_batch.obs['batch'] == batch])
                                for batch in batches}
            batch_max_reads = batches[np.argmax(n_reads_per_cell.values())]

            for batch in batches:
                mult_value = float((np.median(dict_counts_before[batch_max_reads]) / 
                                    np.median(dict_counts_before[batch])))

                adata_batch[adata_batch.obs['batch'] == batch].X *= mult_value

                if 'unspliced' in adata_batch.layers:
                    adata_batch.layers['unspliced'] *= mult_value
                    adata_batch.layers['spliced'] *= mult_value

            dict_counts_after = {batch: adata_batch[adata_batch.obs['batch'] == batch].X.sum(1) for batch in batches}

            fig, axs = plt.subplots(1, 2, figsize=(15,7))
            for batch in batches:
                len_batch = len(dict_counts_before[batch])
                axs[0].plot(np.arange(len_batch)/len_batch, np.log10(np.sort(dict_counts_before[batch])), label=batch)
                axs[1].plot(np.arange(len_batch)/len_batch, np.log10(np.sort(dict_counts_after[batch])))
            axs[0].legend()
            
            
def preprocess_adata_sub(adata_sub, resolution, n_HVG, min_dist=0.8, n_comps=50, seed=10): # MAKE THIS ARGUMENTS COMPULSORY!
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_HVG)
    sc.pl.highly_variable_genes(adata_sub)
    sc.pp.pca(adata_sub, n_comps=n_comps)
    if 'donor' in adata_sub.obs:
        sce.pp.bbknn(adata_sub, metric='angular', batch_key='donor')
    else:
        sc.pp.neighbors(adata_sub, metric='cosine')
    sc.tl.umap(adata_sub, random_state=seed, min_dist=min_dist)
    sc.tl.leiden(adata_sub, random_state=seed, resolution=resolution)