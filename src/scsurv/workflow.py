import torch
import scanpy as sc
import pandas as pd
import numpy as np
from .exp import scSurvExperiment
from .utils import make_inputs, safe_toarray, make_sample_one_hot_mat, vae_results, bulk_deconvolution_results, beta_z_results, spatial_results, optimize_vae, optimize_deepcolor, optimize_scSurv


def run_scSurv(
        sc_adata, bulk_adata, param_save_path, epoch, batch_key, bulk_seed=0, spatial_adata = None, method='efron',
        survival_time_label = 'survival_times', survival_time_censor = 'vital_status', saved_path = None):
    
    x_batch_size_VAE = 1000
    x_batch_size_DeepCOLOR = 1000
    x_batch_size_scSurv = 500
    
    first_lr=0.01
    second_lr=0.01
    third_lr = 0.0001
    bulk_validation_num_or_ratio=0.2
    bulk_test_num_or_ratio=0.2
    
    use_val_loss_mean = True
    model_params = {"z_dim": 20, "h_dim": 100, "num_enc_z_layers": 1, "num_dec_z_layers":1 , "num_dec_p_layers": 1 , "num_dec_b_layers": 1}
    
    patience=10
    usePoisson_sc=True
    x_count, bulk_count = make_inputs(sc_adata, bulk_adata)
    model_params['x_dim'] = x_count.size()[1]
    batch_onehot = make_sample_one_hot_mat(sc_adata, batch_key)

    if saved_path is not None:
        survival_time = None
        cutting_off_0_1 = None
        spatial_count = None
        scsurv_exp = scSurvExperiment(model_params=model_params, x_count=x_count, bulk_count=bulk_count, survival_time=survival_time, cutting_off_0_1=cutting_off_0_1, x_batch_size=x_batch_size_VAE, checkpoint=param_save_path,
                                    usePoisson_sc=usePoisson_sc, batch_onehot=batch_onehot, spatial_count=spatial_count, method=method, use_val_loss_mean=use_val_loss_mean, saved_path=saved_path)
        scsurv_exp.scsurv.load_state_dict(torch.load(saved_path))
        sc_adata.uns['param_save_path'] = saved_path
        model_params_dict = {}
    else:
        survival_time_np = (bulk_adata.obs[survival_time_label]).values
        survival_time = torch.tensor(survival_time_np)

        if survival_time_censor=='vital_status':
            valid_values = ['Dead', 'Alive']
            if not all(value in valid_values for value in bulk_adata.obs[survival_time_censor]):
                raise ValueError("Invalid values found in bulk_adata.obs[survival_time_censor]. Only 'Dead' and 'Alive' are allowed.")
            vital_status_values = np.where(bulk_adata.obs[survival_time_censor] == 'Dead', 1, 0)
            cutting_off_0_1 = torch.tensor(vital_status_values)
        else:
            cutting_off_0_1 = torch.tensor(bulk_adata.obs[survival_time_censor].values)

        model_params_dict = {'1st_lr':first_lr, '2nd_lr':second_lr, '3rd_lr':third_lr, 'patience':patience, 'bulk_seed':bulk_seed,
                            'x_batch_size_VAE':x_batch_size_VAE, 'x_batch_size_DeepCOLOR':x_batch_size_DeepCOLOR, 'x_batch_size_scSurv':x_batch_size_scSurv,
                            'n_var':sc_adata.n_vars, 'usePoisson_sc':usePoisson_sc, 'batch_key':batch_key, 
                            'n_obs_sc':sc_adata.n_obs, 'n_obs_bulk':bulk_adata.n_obs, 'method':method, 'use_val_loss_mean:' : use_val_loss_mean,
                            'bulk_validation_num_or_ratio': bulk_validation_num_or_ratio, 'bulk_test_num_or_ratio': bulk_test_num_or_ratio,}


        if spatial_adata is not None:
            spatial_count = torch.tensor(safe_toarray(spatial_adata.X))
            model_params_dict['n_obs_spatial'] = spatial_adata.n_obs
        else:
            spatial_count = None
        model_params_dict.update(model_params)
        print(model_params_dict)
        scsurv_exp = scSurvExperiment(model_params=model_params, x_count=x_count, bulk_count=bulk_count, survival_time=survival_time, cutting_off_0_1=cutting_off_0_1, x_batch_size=x_batch_size_VAE, checkpoint=param_save_path,
                                    usePoisson_sc=usePoisson_sc, batch_onehot=batch_onehot, spatial_count=spatial_count, method=method, use_val_loss_mean=use_val_loss_mean)
        scsurv_exp = optimize_vae(scsurv_exp=scsurv_exp, first_lr=first_lr, x_batch_size=x_batch_size_VAE, epoch=epoch, patience=patience, param_save_path=param_save_path)
        torch.save(scsurv_exp.scsurv.state_dict(), param_save_path.replace('.pt', '') + '_1st_end.pt')
        scsurv_exp = optimize_deepcolor(scsurv_exp=scsurv_exp, second_lr=second_lr, x_batch_size=x_batch_size_DeepCOLOR, epoch=epoch, patience=patience, param_save_path=param_save_path, spatial_adata=spatial_adata)
        torch.save(scsurv_exp.scsurv.state_dict(), param_save_path.replace('.pt', '') + '_2nd_end.pt')
        scsurv_exp.bulk_data_split(bulk_seed, bulk_validation_num_or_ratio, bulk_test_num_or_ratio, cutting_off_0_1)
        scsurv_exp = optimize_scSurv(scsurv_exp, third_lr = third_lr, x_batch_size=x_batch_size_scSurv, epoch = epoch, patience = patience, param_save_path = param_save_path)    
        torch.save(scsurv_exp.scsurv.state_dict(), param_save_path)
        sc_adata.uns['param_save_path'] = param_save_path
    print('Done scSurv')
    return sc_adata, bulk_adata, model_params_dict, spatial_adata, scsurv_exp

def post_process(scsurv_exp, sc_adata, bulk_adata, spatial_adata=None, save_memory=False):
    if save_memory==False:
        sc_adata = vae_results(scsurv_exp, sc_adata)
    sc_adata, bulk_adata = bulk_deconvolution_results(scsurv_exp, sc_adata, bulk_adata, save_memory)
    sc_adata = beta_z_results(scsurv_exp, sc_adata)
    if spatial_adata is not None:
        sc_adata, spatial_adata = spatial_results(scsurv_exp, sc_adata, spatial_adata)
    print('Done post process')
    return sc_adata, bulk_adata, spatial_adata

def scSurv_preprocess(sc_adata, bulk_adata, per =  0.01, n_top_genes = 5000, highly_variable='bulk'):
        
        common_genes_before = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
        sc_adata = sc_adata[:, common_genes_before].copy()
        bulk_adata = bulk_adata[:, common_genes_before].copy()
        print('common_genes_before', len(common_genes_before))

        sc_min_cells = int(sc_adata.n_obs * per)
        bulk_min_cells = int(bulk_adata.n_obs * per)
        sc.pp.filter_genes(sc_adata, min_cells=sc_min_cells)
        sc.pp.filter_genes(bulk_adata, min_cells=bulk_min_cells)

        common_genes_filtered = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
        sc_adata = sc_adata[:, common_genes_filtered].copy()
        bulk_adata = bulk_adata[:, common_genes_filtered].copy()
        print('common_genes_filtered', len(common_genes_filtered))

        raw_sc_adata = sc_adata.copy()
        raw_bulk_adata = bulk_adata.copy()
        if highly_variable=='bulk':
            sc.pp.normalize_total(bulk_adata)
            sc.pp.log1p(bulk_adata)
            sc.pp.highly_variable_genes(bulk_adata, n_top_genes=n_top_genes)
            bulk_adata = bulk_adata[:, bulk_adata.var['highly_variable']].copy()
        elif highly_variable=='sc':
            sc.pp.normalize_total(sc_adata)
            sc.pp.log1p(sc_adata)
            sc.pp.highly_variable_genes(sc_adata, n_top_genes=n_top_genes)
            sc_adata = sc_adata[:, sc_adata.var['highly_variable']].copy()

        common_genes_highly_variable = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
        if len(common_genes_highly_variable) == 0:
            raise ValueError("No common genes found between sc_adata and bulk_adata.")
        sc_adata = sc_adata[:, common_genes_highly_variable].copy()
        bulk_adata = bulk_adata[:, common_genes_highly_variable].copy()
        print('common_genes_highly_variable', len(common_genes_highly_variable))
        
        sc_adata.X = raw_sc_adata[:, sc_adata.var_names].X
        bulk_adata.X = raw_bulk_adata[:, bulk_adata.var_names].X
        return sc_adata, bulk_adata

def scSurv_preprocess_spatial(sc_adata, bulk_adata, spatial_adata, per =  0.01, n_top_genes = 5000, highly_variable='bulk'):
        common_genes_before = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
        common_genes_before = np.intersect1d(common_genes_before, spatial_adata.var_names)
        sc_adata = sc_adata[:, common_genes_before].copy()
        bulk_adata = bulk_adata[:, common_genes_before].copy()
        spatial_adata = spatial_adata[:, common_genes_before].copy()
        print('common_genes_before', len(common_genes_before))

        sc_min_cells = int(sc_adata.n_obs * per)
        bulk_min_cells = int(bulk_adata.n_obs * per)
        sp_min_cells = int(spatial_adata.n_obs * per)
        sc.pp.filter_genes(sc_adata, min_cells=sc_min_cells)
        sc.pp.filter_genes(bulk_adata, min_cells=bulk_min_cells)
        sc.pp.filter_genes(spatial_adata, min_cells=sp_min_cells)
        common_genes_filtered = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
        common_genes_filtered = np.intersect1d(common_genes_filtered, spatial_adata.var_names)
        sc_adata = sc_adata[:, common_genes_filtered].copy()
        bulk_adata = bulk_adata[:, common_genes_filtered].copy()
        spatial_adata = spatial_adata[:, common_genes_filtered].copy()
        print('common_genes_filtered', len(common_genes_filtered))

        raw_sc_adata = sc_adata.copy()
        raw_bulk_adata = bulk_adata.copy()
        raw_spatial_adata = spatial_adata.copy()
        if highly_variable=='bulk':
            sc.pp.normalize_total(bulk_adata)
            sc.pp.log1p(bulk_adata)
            sc.pp.highly_variable_genes(bulk_adata, n_top_genes=n_top_genes)
            bulk_adata = bulk_adata[:, bulk_adata.var['highly_variable']].copy()
        elif highly_variable=='spatial':
            sc.pp.normalize_total(spatial_adata)
            sc.pp.log1p(spatial_adata)
            sc.pp.highly_variable_genes(spatial_adata, n_top_genes=n_top_genes)
            spatial_adata = spatial_adata[:, spatial_adata.var['highly_variable']].copy()
        elif highly_variable=='sc':
            sc.pp.normalize_total(sc_adata)
            sc.pp.log1p(sc_adata)
            sc.pp.highly_variable_genes(sc_adata, n_top_genes=n_top_genes)
            sc_adata = sc_adata[:, sc_adata.var['highly_variable']].copy()

        common_genes_highly_variable = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
        common_genes_highly_variable = np.intersect1d(common_genes_highly_variable, spatial_adata.var_names)
        if len(common_genes_highly_variable) == 0:
            raise ValueError("No common genes found between sc_adata, bulk_adata and spatial_adata.")
        sc_adata = sc_adata[:, common_genes_highly_variable].copy()
        bulk_adata = bulk_adata[:, common_genes_highly_variable].copy()
        spatial_adata = spatial_adata[:, common_genes_highly_variable].copy()
        print('common_genes_highly_variable', len(common_genes_highly_variable))

        sc_adata.X = raw_sc_adata[:, sc_adata.var_names].X
        bulk_adata.X = raw_bulk_adata[:, bulk_adata.var_names].X
        spatial_adata.X = raw_spatial_adata[:, spatial_adata.var_names].X
        return sc_adata, bulk_adata, spatial_adata
