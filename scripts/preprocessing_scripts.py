import numpy as np
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
import scrublet as scr
import scipy.sparse as spr


def cell_selection_plotting(adata, cell_count_options, palette, cmap):
    dict_argidx = {}
    
    adata_cr_counts = adata.copy()
    
    for cell_counts in cell_count_options:
        n_cells = min(cell_counts, int(0.8 * adata_cr_counts.X.shape[0]))
        n_counts = np.asarray(adata_cr_counts.X.sum(1)).flatten()
        cutoff = np.sort(n_counts.copy())[::-1][n_cells]
        argidx = np.argwhere(n_counts >= cutoff).flatten()
        dict_argidx[cell_counts] = adata_cr_counts.obs.index[argidx]

    adata_cr_counts = adata_cr_counts[dict_argidx[cell_count_options[0]],:] # select the biggest number

    sc.pp.filter_genes(adata_cr_counts, min_cells=2)
    sc.pp.log1p(adata_cr_counts)
    adata_cr_counts.obs['n_counts'] = np.log10(adata_cr_counts.X.sum(1))
    sc.tl.pca(adata_cr_counts)
    sc.pp.neighbors(adata_cr_counts)
    sc.tl.umap(adata_cr_counts)

    adata_cr_counts.obs['n_cells'] = cell_count_options[-1] # select the fewest cells

    for idx in range(len(cell_count_options[:-1])):
        current_argidx, next_argidx = dict_argidx[cell_count_options[idx]], dict_argidx[cell_count_options[idx + 1]]
        idx_cells = current_argidx[~np.isin(current_argidx, next_argidx)]
        adata_cr_counts.obs['n_cells'].loc[idx_cells] = cell_count_options[idx]

    adata_cr_counts.obs['n_cells'] = adata_cr_counts.obs['n_cells'].astype(str)  
    
    sc.pl.umap(adata_cr_counts, color=['n_counts', 'n_cells', 'MALAT1'], palette=palette, cmap=cmap, ncols=3, alpha=0.7)
    sns.jointplot(np.asarray(adata_cr_counts[:, 'MALAT1'].X.todense()).flatten(), adata_cr_counts.obs['n_counts'].values)

    return adata_cr_counts


def cell_selection(adata, selected_n_cells, MALAT1_threshold):
    n_counts = np.asarray(adata.X.sum(1)).flatten()
    cutoff = np.sort(n_counts.copy())[::-1][selected_n_cells]
    argidx_ncells = np.where(n_counts >= cutoff)[0]
    argidx_MALAT1 = np.argwhere(adata[:, 'MALAT1'].X > MALAT1_threshold)[:,0]
    argidx_both = np.intersect1d(argidx_MALAT1, argidx_ncells)
    adata = adata[argidx_both]
    n_cells = len(argidx_both)
    
    sc.pp.filter_genes(adata, min_cells=3)
    
    print('Number of median counts per cell: {}'.format(np.median(np.asarray(adata.X.sum(1)[:,None]))))
    
    return adata


def plot_mito(adata):
    mito_genes = adata.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    sc.pp.calculate_qc_metrics(adata, inplace = True)
    
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', ax=axs[0], show=False)
    sc.pl.scatter(adata, x='n_counts', y='n_genes_by_counts', ax=axs[1])
    
    return adata


def run_scrublet(adata, min_dist):
    scrub = scr.Scrublet(spr.csc_matrix(adata.X))
    doublet_scores, predicted_doublets = scrub.scrub_doublets(distance_metric='cosine', n_prin_comps=15,)

    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=min_dist))
    scrub.plot_embedding('UMAP', order_points=True);
    
    adata.obs['scrublet_doublet'] = predicted_doublets
    
    return adata


def prepare_norm(adata, n_comps=50, resolution=0.1, counts_per_cell_after=1e6, seed=10):
    # We create a mock dataset to obtain the vector for Z types
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=counts_per_cell_after)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=n_comps)
    sc.pp.neighbors(adata_pp, metric='cosine')
    sc.tl.umap(adata_pp, random_state=seed)
    sc.tl.leiden(adata_pp, key_added='groups', resolution=resolution, random_state=seed)
    sc.pl.umap(adata_pp, color='groups')

    #Preprocess variables for scran normalization
    input_groups = adata_pp.obs['groups'].values.astype(int)
    data_mat = adata.X.todense().T 
    
    return data_mat, input_groups


def apply_norm(adata, size_factors):
    adata.raw = adata
    adata.obs['size_factors'] = size_factors
    
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    sc.pl.scatter(adata, 'size_factors', 'n_counts', ax=axs[0], show=False)
    #Keep the count data in a counts layer
    adata.layers["raw_counts"] = adata.X.copy()
    
    sns.distplot(size_factors, bins=50, kde=False, ax=axs[1])
    
    adata = adata[adata.obs['scrublet_doublet'] == False]
    
    #Normalize adata 
    adata.X /= adata.obs['size_factors'].values[:,None]
    
    return adata


def adata_preprocessing(adata, sample, hvg_min_mean, hvg_max_mean, hvg_min_disp, seed, leiden_resolution, dir_adata_save, cmap, palette):
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, metric='cosine', n_neighbors=int(0.5 * len(adata) ** 0.5))
    sc.tl.umap(adata, random_state = seed)
    sc.tl.leiden(adata, resolution=leiden_resolution, random_state = seed)
    
    sc.pl.umap(adata, color=['leiden', 'log1p_total_counts'], 
           cmap=cmap, palette=palette)
    
    adata.write_h5ad(f'{dir_adata_save}/{sample}_processed.h5')