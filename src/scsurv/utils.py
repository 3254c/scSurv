import torch
import pandas as pd
import numpy as np

def safe_toarray(x):
    if type(x) != np.ndarray:
        x = x.toarray()
        if not np.all(x == np.floor(x)):
            raise ValueError('target layer of adata should be raw count')
        return x
    else:
        if not np.all(x == np.floor(x)):
            raise ValueError('target layer of adata should be raw count')
        return x

def make_sample_one_hot_mat(adata, sample_key):
    print('make_sample_one_hot_mat')
    if sample_key is not None:
        sidxs = np.sort(adata.obs[sample_key].unique())
        b = np.array([
            (sidxs == sidx).astype(int)
            for sidx in adata.obs[sample_key]]).astype(float)
        b = torch.tensor(b).float()
    else:
        # 観測値の数に基づいてゼロ配列を作成（形状が(観測値の数, 1)）
        b = np.zeros((len(adata.obs_names), 1))
        b = torch.tensor(b).float()
    return b

def input_checks(adata, layer_name):
    if layer_name == 'X':
        if np.sum((adata.X - adata.X.astype(int)))**2 != 0:
            raise ValueError('`X` includes non integer number, while count data is required for `X`.')
    else:
        if np.sum((adata.layers[layer_name] - adata.layers[layer_name].astype(int)))**2 != 0:
            raise ValueError(f'layers `{layer_name}` includes non integer number, while count data is required for `{layer_name}`.')

def make_inputs(sc_adata, bulk_adata, layer_name='X'):
    input_checks(sc_adata, layer_name)
    if layer_name == 'X':
        x = torch.tensor(safe_toarray(sc_adata.X))
        if bulk_adata is not None:
            s = torch.tensor(safe_toarray(bulk_adata.X))
        else:
            s = None
    else:
        x = torch.tensor(safe_toarray(sc_adata.layers[layer_name]))
        if bulk_adata is not None:
            s = torch.tensor(safe_toarray(bulk_adata.layers[layer_name]))
        else:
            s = None
    return x, s

def optimize_vae(scsurv_exp, first_lr, x_batch_size, epoch, patience, param_save_path):
    print('Start first opt', 'lr=', first_lr)
    scsurv_exp.scsurv.sc_mode()
    scsurv_exp.initialize_optimizer(first_lr)
    scsurv_exp.initialize_loader(x_batch_size)
    stop_epoch_vae = scsurv_exp.train_total(epoch, patience)
    scsurv_exp.scsurv.load_state_dict(torch.load(param_save_path))
    val_loss_vae = scsurv_exp.evaluate(mode='val')
    test_loss_vae = scsurv_exp.evaluate(mode='test')    
    print(f'Done {scsurv_exp.scsurv.mode} mode,', f'Val Loss: {val_loss_vae}', f'Test Loss: {test_loss_vae}')
    return scsurv_exp

def optimize_deepcolor(scsurv_exp, second_lr, x_batch_size, epoch, patience, param_save_path, spatial_adata):
    scsurv_exp.scsurv.bulk_mode()
    scsurv_exp.initialize_optimizer(second_lr)
    scsurv_exp.initialize_loader(x_batch_size)
    print(f'{scsurv_exp.scsurv.mode} mode', 'lr=', second_lr)
    stop_epoch_bulk = scsurv_exp.train_total(epoch, patience)
    scsurv_exp.scsurv.load_state_dict(torch.load(param_save_path))
    val_loss_bulk = scsurv_exp.evaluate(mode='val')
    test_loss_bulk = scsurv_exp.evaluate(mode='test')
    print(f'Done {scsurv_exp.scsurv.mode} mode,', f'Val Loss: {val_loss_bulk}', f'Test Loss: {test_loss_bulk}')
    if spatial_adata is not None:
        scsurv_exp.scsurv.spatial_mode()
        scsurv_exp.initialize_optimizer(second_lr)
        scsurv_exp.initialize_loader(x_batch_size)
        print(f'{scsurv_exp.scsurv.mode} mode', 'lr=', second_lr)
        stop_epoch_spatial = scsurv_exp.train_total(epoch, patience)
        scsurv_exp.scsurv.load_state_dict(torch.load(param_save_path))
        val_loss_spatial = scsurv_exp.evaluate(mode='val')
        test_loss_spatial = scsurv_exp.evaluate(mode='test')
        print(f'Done {scsurv_exp.scsurv.mode} mode,', f'Val Loss: {val_loss_spatial}', f'Test Loss: {test_loss_spatial}')
    return scsurv_exp

def optimize_scSurv(scsurv_exp, third_lr, x_batch_size, epoch, patience, param_save_path):
    scsurv_exp.scsurv.hazard_beta_z_mode()
    scsurv_exp.initialize_optimizer(third_lr)
    scsurv_exp.initialize_loader(x_batch_size)
    print(f'{scsurv_exp.scsurv.mode} mode', 'lr=', third_lr)
    stop_epoch_beta_z = scsurv_exp.train_total(epoch, patience)
    scsurv_exp.scsurv.load_state_dict(torch.load(param_save_path))
    train_loss_beta_z = abs(scsurv_exp.evaluate_train())
    val_loss_beta_z = abs(scsurv_exp.evaluate(mode='val'))
    test_loss_beta_z = abs(scsurv_exp.evaluate(mode='test'))
    print(f'Done {scsurv_exp.scsurv.mode} mode,', f'Train c-index: {train_loss_beta_z}', f'Val c-index: {val_loss_beta_z}', f'Test c-index: {test_loss_beta_z}')
    return scsurv_exp


def vae_results(scsurv_exp, sc_adata):
    print('vae_results')
    with torch.no_grad():
        torch.cuda.empty_cache()
        batch_onehot = scsurv_exp.x_data_manager.batch_onehot.to(scsurv_exp.device)
        batch_onehot_np = batch_onehot.detach().cpu().numpy()
        sc_adata.obsm['batch_onehot'] = batch_onehot_np
        x = scsurv_exp.x_data_manager.x_count.to(scsurv_exp.device)
        x_np = x.detach().cpu().numpy()
        xb = torch.cat([x, batch_onehot], dim=-1)
        z, qz = scsurv_exp.scsurv.enc_z(xb)
        zl = qz.loc
        xld_mean = None
        for i in range(100):
            zzz = qz.sample()
            zb = torch.cat([zzz, batch_onehot], dim=-1)
            xxx_np = scsurv_exp.scsurv.dec_z2x(zb).detach().cpu().numpy()
            if xld_mean is None:
                xld_mean = xxx_np
            else:
                xld_mean += (xxx_np - xld_mean) / (i + 1)
            del zzz, zb, xxx_np
            torch.cuda.empty_cache()
        xld_np = xld_mean
        sc_adata.obsm['zl'] = zl.detach().cpu().numpy()
        sc_adata.layers['xld'] = xld_np
        xnorm_mat=scsurv_exp.x_data_manager.xnorm_mat
        xnorm_mat_np = xnorm_mat.cpu().detach().numpy()
        x_df = pd.DataFrame(x_np, columns=list(sc_adata.var_names))
        xld_df = pd.DataFrame(xld_np,columns=list(sc_adata.var_names))
        train_idx = scsurv_exp.x_data_manager.train_idx
        val_idx = scsurv_exp.x_data_manager.validation_idx
        test_idx = scsurv_exp.x_data_manager.test_idx

        x_correlation_gene=(xld_df).corrwith(x_df / xnorm_mat_np).mean()
        train_x_correlation_gene = (xld_df.T[train_idx].T).corrwith((x_df / xnorm_mat_np).T[train_idx].T).mean()
        val_x_correlation_gene = (xld_df.T[val_idx].T).corrwith((x_df / xnorm_mat_np).T[val_idx].T).mean()
        test_x_correlation_gene = (xld_df.T[test_idx].T).corrwith((x_df / xnorm_mat_np).T[test_idx].T).mean()

        print('all_x_correlation_gene', f"{x_correlation_gene:.3f}", 'train_x_correlation_gene', f"{train_x_correlation_gene:.3f}", 'val_x_correlation_gene', f"{val_x_correlation_gene:.3f}", 'test_x_correlation_gene', f"{test_x_correlation_gene:.3f}")
        return sc_adata

def bulk_deconvolution_results(scsurv_exp, sc_adata, bulk_adata, save_memory):
    print('deconvolution_results')
    with torch.no_grad():
        torch.cuda.empty_cache()
        batch_onehot = scsurv_exp.x_data_manager.batch_onehot.to(scsurv_exp.device)
        x = scsurv_exp.x_data_manager.x_count.to(scsurv_exp.device)
        xb = torch.cat([x, batch_onehot], dim=-1)
        z, qz = scsurv_exp.scsurv.enc_z(xb)
        p_mean = None
        for i in range(100):
            zzz = qz.sample()
            ppp = scsurv_exp.scsurv.dec_z2p_bulk(zzz).detach().cpu().numpy()
            if p_mean is None:
                p_mean = ppp
            else:
                p_mean += (ppp - p_mean) / (i + 1)
            del zzz, ppp
            torch.cuda.empty_cache()
        bulk_pl_np = p_mean
        sc_adata.obsm['map2bulk'] = bulk_pl_np.transpose()
        if save_memory==False:
            bulk_scoeff_np = scsurv_exp.scsurv.softplus(scsurv_exp.scsurv.log_bulk_coeff).cpu().detach().numpy()
            bulk_scoeff_add_np = scsurv_exp.scsurv.softplus(scsurv_exp.scsurv.log_bulk_coeff_add).cpu().detach().numpy()
            xld_np = sc_adata.layers['xld']
            bulk_hat_np = np.matmul(bulk_pl_np, xld_np * bulk_scoeff_np) + bulk_scoeff_add_np #spatialの発現のmean parameter #matmulは行列の積 #.cpu().detach().numpy()
            # bulk_adata.obsm['map2sc'] = bulk_p_df.transpose().values
            if bulk_adata is not None:
                bulk_hat_df = pd.DataFrame(bulk_hat_np, columns=list(sc_adata.var_names))
                bulk_norm_mat=scsurv_exp.bulk_data_manager.bulk_norm_mat
                bulk_norm_mat_np = bulk_norm_mat.cpu().detach().numpy()
                bulk_count = scsurv_exp.bulk_data_manager.bulk_count
                bulk_count_df = pd.DataFrame(bulk_count.cpu().detach().numpy(), columns=list(sc_adata.var_names))
                bulk_adata.layers['bulk_hat'] = pd.DataFrame(bulk_hat_np, index = list(bulk_adata.obs_names), columns=list(sc_adata.var_names))
                bulk_correlation_gene=(bulk_hat_df).corrwith(bulk_count_df / bulk_norm_mat_np).mean()
                print('bulk_correlation_gene', bulk_correlation_gene)
        return sc_adata, bulk_adata

def spatial_results(scsurv_exp, sc_adata, spatial_adata):
    print('spatial_results')
    with torch.no_grad():
        torch.cuda.empty_cache()
        batch_onehot = scsurv_exp.x_data_manager.batch_onehot.to(scsurv_exp.device)
        x = scsurv_exp.x_data_manager.x_count.to(scsurv_exp.device)
        xb = torch.cat([x, batch_onehot], dim=-1)
        z, qz = scsurv_exp.scsurv.enc_z(xb)
        ppp_list = []
        for _ in range(100):
            zzz = qz.sample()
            ppp = scsurv_exp.scsurv.dec_z2p_spatial(zzz).detach().cpu().numpy()
            ppp_list.append(ppp)
        spatial_pl_np = np.mean(ppp_list, axis=0) #scRNAseq中の各細胞がspatialやbulkに占める割合
        del zzz, ppp
        #spatial_pl_np = scsurv_exp.scsurv.dec_z2p_spatial(zl).detach().cpu().numpy()
        spatial_coeff_np = scsurv_exp.scsurv.softplus(scsurv_exp.scsurv.log_spatial_coeff).cpu().detach().numpy()
        spatial_coeff_add_np = scsurv_exp.scsurv.softplus(scsurv_exp.scsurv.log_spatial_coeff_add).cpu().detach().numpy()
        xld_np = sc_adata.layers['xld']
        spatial_hat_np = np.matmul(spatial_pl_np, xld_np * spatial_coeff_np) + spatial_coeff_add_np #spatialの発現のmean parameter #matmulは行列の積 #.cpu().detach().numpy()
        spatial_p_df = pd.DataFrame(spatial_pl_np.transpose(), index=sc_adata.obs_names, columns=spatial_adata.obs_names)
        sc_adata.obsm['map2spatial'] = spatial_p_df.values
        # spatial_adata.obsm['map2sc'] = spatial_p_df.transpose().values
        beta_z_p_spatial = sc_adata.obsm['map2spatial'] * sc_adata.obs['beta_z'].values[:, None]
        spatial_adata.obs['exp_beta_p_spatial_sum'] = np.exp(beta_z_p_spatial.sum(axis=0))
        spatial_adata.obs['Hazard_rates'] = spatial_adata.obs['exp_beta_p_spatial_sum'] / spatial_adata.obs['exp_beta_p_spatial_sum'].mean()
        spatial_norm_mat=scsurv_exp.spatial_data_manager.spatial_norm_mat
        spatial_norm_mat_np = spatial_norm_mat.cpu().detach().numpy()
        spatial_count = scsurv_exp.spatial_data_manager.spatial_count
        spatial_count_df = pd.DataFrame(spatial_count.cpu().detach().numpy(), columns=list(spatial_adata.var_names))
        spatial_hat_df = pd.DataFrame(spatial_hat_np,columns=list(spatial_adata.var_names))
        spatial_adata.layers['spatial_hat'] = pd.DataFrame(spatial_hat_np, index = list(spatial_adata.obs_names), columns=list(spatial_adata.var_names))
        spatial_correlation_gene=(spatial_hat_df).corrwith(spatial_count_df / spatial_norm_mat_np).mean()
        print('spatial_correlation_gene', spatial_correlation_gene)
        return sc_adata, spatial_adata

def beta_z_results(scsurv_exp, sc_adata):
    with torch.no_grad():
        torch.cuda.empty_cache()
        batch_onehot = scsurv_exp.x_data_manager.batch_onehot.to(scsurv_exp.device)
        x = scsurv_exp.x_data_manager.x_count.to(scsurv_exp.device)
        xb = torch.cat([x, batch_onehot], dim=-1)
        z, qz = scsurv_exp.scsurv.enc_z(xb)
        beta_list = []
        zl = qz.loc
        for _ in range(100):
            z_sample = qz.sample()
            b_sample= scsurv_exp.scsurv.dec_beta_z(z_sample).detach().cpu().numpy()
            beta_list.append(b_sample)
        beta_z_np = np.mean(beta_list, axis=0) #scRNAseq中の各細胞がspatialやbulkに占める割合
        beta_zl_np = scsurv_exp.scsurv.dec_beta_z(zl).detach().cpu().numpy()
        sc_adata.obs['raw_beta_z'] = beta_z_np
        sc_adata.obs['raw_beta_zl'] = beta_zl_np
        sc_adata.obs['beta_z'] = beta_z_np - np.min(beta_z_np)
        return sc_adata