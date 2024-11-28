import src.scsurv.workflow as workflow
import pandas as pd
import scanpy as sc

sc_adata = sc.read_h5ad('/home/mizukoshi/scSurv/simulated_BRCA_5000genes_CID4471_celltype_minor_300_major_reso0.1_sc_reference.h5ad')
bulk_adata = sc.read_h5ad(f'/home/mizukoshi/scSurv/single_beta_bulk_BRCA_betascale5_alive0_major_reso0.1_0.h5ad')
bulk_ncells_per_leiden = bulk_adata.uns['n_cells_per_celltype_df'].T
setting_beta = pd.Series(bulk_adata.uns['setting_beta'], index = bulk_ncells_per_leiden.columns)
cell_type_df = sc_adata.obs['celltype_minor']
sc_adata.obs['setting_beta'] = sc_adata.obs['celltype_minor'].map(setting_beta).astype(float)

batch_key = 'orig.ident'
exp_name = f'tutolial_simulation.pt'
epoch = 10000

sc_adata, bulk_adata, model_params_dict, spatial_adata, vaesm_exp = workflow.run_scSurv(sc_adata, bulk_adata, exp_name, epoch, batch_key)
sc.pl.umap(sc_adata, color = ['celltype_minor'], save = '_celltype_minor.png')
sc.pl.umap(sc_adata, color = ['setting_beta'], save = '_setting_beta.png')
sc.pl.umap(sc_adata, color = ['beta_z'], save = '_estimated_beta.png')

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
print('scSurv beta correlation', beta_corr)