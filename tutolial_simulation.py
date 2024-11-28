#The scRNA-seq data used in this simulation is derived from the study by Wu et al. (https://doi.org/10.1038/s41588-021-00911-1).
import src.scsurv.workflow as workflow
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sc_adata = sc.read_h5ad('/home/mizukoshi/scSurv/tutolial_simulated_sc_reference.h5ad')
bulk_adata = sc.read_h5ad('/home/mizukoshi/scSurv/tutolial_simulated_bulk.h5ad')
bulk_ncells_per_leiden = bulk_adata.uns['n_cells_per_celltype_df'].T
setting_beta = pd.Series(bulk_adata.uns['setting_beta'], index = bulk_ncells_per_leiden.columns)
cell_type_df = sc_adata.obs['celltype_minor']
sc_adata.obs['setting_beta'] = sc_adata.obs['celltype_minor'].map(setting_beta).astype(float)

batch_key = 'orig.ident'
exp_name = f'tutolial_simulation.pt'
epoch = 10000

sc_adata, bulk_adata, model_params_dict, spatial_adata, vaesm_exp = workflow.run_scSurv(sc_adata, bulk_adata, exp_name, epoch, batch_key)

fig, ax = plt.subplots(figsize=(14, 6))
sc.pl.umap(sc_adata, color='celltype_minor', ax=ax, show=False, colorbar_loc = None, frameon = False, title='celltype_minor')
fig.subplots_adjust(right=0.5)
plt.savefig('umap_celltype_minor.png')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(sc_adata, color='setting_beta', ax=ax1, show=False, colorbar_loc=None, frameon=False, title='Setting contribution')
sc.pl.umap(sc_adata, color='beta_z', ax=ax2, show=False, colorbar_loc=None, frameon=False, title='Estimated contribution')
plt.savefig('umap_setting_and_estimated_beta.png')

estimated_p_sc_level = pd.DataFrame(sc_adata.obsm['map2bulk'].T, index=bulk_adata.obs_names, columns=sc_adata.obs_names)
bulk_ncells = bulk_adata.uns['n_cells_per_celltype_df'].T
estimated_p_cluster = pd.DataFrame(index=bulk_ncells.index, columns=bulk_ncells.columns)
for item in bulk_ncells.columns:
    item_idx = sc_adata[sc_adata.obs['celltype_minor']==item].obs.index
    estimated_p_cluster[item] = estimated_p_sc_level[item_idx].sum(axis=1)

deconv_corr = estimated_p_cluster.corrwith(bulk_ncells, axis=1, method='pearson')
print('scSurv deconvolution correlation', deconv_corr.median())

beta_means = {item: sc_adata[sc_adata.obs['celltype_minor'] == item].obs['beta_z'].mean() 
             for item in bulk_ncells.columns}
beta_corr = pd.Series(beta_means).corr(setting_beta)
print('scSurv contribution correlation', beta_corr)