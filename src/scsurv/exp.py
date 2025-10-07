import torch
from .modules import scSurv
from .dataset import ScDataManager, BulkDataManager, SpatialDataManager
from statistics import mean
import torch.distributions as dist
import math
from collections import deque
from lifelines.utils import concordance_index


class EarlyStopping:
    def __init__(self, patience, path):
        self.patience = patience
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.path = path

    def __call__(self, val_loss, model):
        if self.best_score is None:
            self.best_score = val_loss
            self.checkpoint(model)
        elif val_loss > self.best_score:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = val_loss
            self.checkpoint(model)
            self.counter = 0
            
    def checkpoint(self, model):
        torch.save(model.state_dict(), self.path)

class scSurvExperiment:
    def __init__(self, model_params, x_count, bulk_count, survival_time, cutting_off_0_1, x_batch_size, checkpoint, usePoisson_sc, 
                batch_onehot, spatial_count, method, use_val_loss_mean, saved_path=None):
        print('torch.cuda.is_available()', torch.cuda.is_available())
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.batch_onehot = batch_onehot
        self.device = torch.device(device)
        self.x_data_manager = ScDataManager(x_count, batch_size=x_batch_size, batch_onehot=batch_onehot)
        self.spatial_count = spatial_count
        self.model_params = model_params
        self.checkpoint=checkpoint
        self.usePoisson_sc = usePoisson_sc
        self.epoch = 0
        self.method = method
        self.use_val_loss_mean = use_val_loss_mean
        self.bulk_test_num_or_ratio = None
        self.bulk_validation_num_or_ratio = None

        if saved_path is not None:
            saved_param = torch.load(saved_path, map_location='cpu')
            bulk_num = saved_param['dec_z2p_bulk.h2ld.weight'].shape[0]
            spatial_num = 0
            self.scsurv = scSurv(bulk_num = bulk_num, spatial_num = spatial_num, batch_onehot_dim = batch_onehot.shape[1], **self.model_params)
        else:
            self.bulk_data_manager = BulkDataManager(bulk_count, survival_time = survival_time, cutting_off_0_1=cutting_off_0_1)
            self.bulk_count = self.bulk_data_manager.bulk_count.to(self.device)
            self.bulk_norm_mat = self.bulk_data_manager.bulk_norm_mat.to(self.device)        
            self.cutting_off_0_1 = self.bulk_data_manager.cutting_off_0_1.to(self.device)
            self.survival_time = self.bulk_data_manager.survival_time.to(self.device)
            if spatial_count is not None:
                self.spatial_data_manager = SpatialDataManager(spatial_count)
                self.spatial_count = self.spatial_data_manager.spatial_count.to(self.device)
                self.spatial_norm_mat = self.spatial_data_manager.spatial_norm_mat.to(self.device)
                spatial_num = self.spatial_data_manager.spatial_count.shape[0]    
            else:
                self.spatial_norm_mat = None
                spatial_num = 0
            self.scsurv = scSurv(bulk_num = self.bulk_data_manager.bulk_count.shape[0], spatial_num = spatial_num, batch_onehot_dim = batch_onehot.shape[1], **self.model_params)
        self.scsurv.to(self.device)


    def bulk_data_split(self, n_bulk_split, validation_num, test_num, censor_np):
        self.bulk_test_num_or_ratio = test_num
        self.bulk_validation_num_or_ratio = validation_num
        self.bulk_data_manager.bulk_split(n_bulk_split, validation_num, test_num, censor_np)

    def elbo_loss(self, x, xnorm_mat, bulk_count, bulk_norm_mat, spatial_count, spatial_norm_mat, batch_onehot, bulk_idx):
        z, qz, x_hat, p_bulk, p_spatial, bulk_hat, spatial_hat, theta_x, theta_bulk, theta_spatial, beta_z = self.scsurv(x, batch_onehot)
        if self.scsurv.mode == 'sc':
            elbo_loss = self.calc_kld(qz).sum()
            if self.usePoisson_sc:
                elbo_loss += self.calc_poisson_loss(ld=x_hat, norm_mat=xnorm_mat, obs=x).sum()
            else:
                elbo_loss += self.calc_nb_loss(x_hat, xnorm_mat, theta_x, x).sum()
        elif self.scsurv.mode == 'bulk':
            elbo_loss = self.calc_nb_loss(bulk_hat, bulk_norm_mat, theta_bulk, bulk_count).sum()
        elif self.scsurv.mode == 'spatial':
            elbo_loss = self.calc_nb_loss(spatial_hat, spatial_norm_mat, theta_spatial, spatial_count).sum()
        elif self.scsurv.mode == 'beta_z' or self.scsurv.mode == 'gamma_z':
            time_indicators = self.cutting_off_0_1.to(self.device)
            survival_time = self.survival_time.to(self.device)
            linear_predictor = torch.sum(beta_z * p_bulk, dim=1)
            neg_log_partial_likelihood, sorted_survival_time, sorted_linear_predictor, sorted_event_observed = self.calc_hazard_loss_beta_z(linear_predictor, survival_time, time_indicators, bulk_idx, self.method)
            elbo_loss = neg_log_partial_likelihood
            if not self.scsurv.training:
                c_index = concordance_index(sorted_survival_time.detach().cpu().numpy(), -sorted_linear_predictor.detach().cpu().numpy(), sorted_event_observed.detach().cpu().numpy())
                elbo_loss = - c_index * x.shape[0]
        return(elbo_loss)
        
    def train_epoch(self):
        self.scsurv.train()
        total_loss = 0
        entry_num = 0
        for x, xnorm_mat, batch_onehot in self.x_data_manager.train_loader:
            x = x.to(self.device)
            xnorm_mat = xnorm_mat.to(self.device)
            batch_onehot = batch_onehot.to(self.device)
            self.scsurv_optimizer.zero_grad()
            bulk_idx = self.bulk_data_manager.train_idx.to(self.device)
            loss = self.elbo_loss(x, xnorm_mat, self.bulk_count, self.bulk_norm_mat,  self.spatial_count, self.spatial_norm_mat, batch_onehot, bulk_idx)
            loss.backward()
            self.scsurv_optimizer.step()
            entry_num += x.shape[0]
            total_loss = total_loss + loss.item()
        loss_val = total_loss / entry_num
        return(loss_val)
        
    def evaluate(self, mode='test'):
        with torch.no_grad():
            self.scsurv.eval()
            if mode == 'test':            
                x = self.x_data_manager.test_x.to(self.device)
                xnorm_mat = self.x_data_manager.test_xnorm_mat.to(self.device)
                batch_onehot = self.x_data_manager.test_batch_onehot.to(self.device)
                bulk_idx = self.bulk_data_manager.test_idx.to(self.device)
            else:
                x = self.x_data_manager.validation_x.to(self.device)
                xnorm_mat = self.x_data_manager.validation_xnorm_mat.to(self.device)
                batch_onehot = self.x_data_manager.validation_batch_onehot.to(self.device)
                bulk_idx = self.bulk_data_manager.validation_idx.to(self.device)
            loss = self.elbo_loss(x, xnorm_mat, self.bulk_count, self.bulk_norm_mat,  self.spatial_count, self.spatial_norm_mat, batch_onehot, bulk_idx)
            entry_num = x.shape[0]
            loss_val = loss / entry_num
            return(loss_val.item())
    
    def evaluate_train(self):
        with torch.no_grad():
            self.scsurv.eval()
            x = self.x_data_manager.train_x.to(self.device)
            xnorm_mat = self.x_data_manager.train_xnorm_mat.to(self.device)
            batch_onehot = self.x_data_manager.train_batch_onehot.to(self.device)
            idx = self.bulk_data_manager.train_idx.to(self.device)
            loss = self.elbo_loss(x, xnorm_mat, self.bulk_count, self.bulk_norm_mat,  self.spatial_count, self.spatial_norm_mat, batch_onehot, idx)
            entry_num = x.shape[0]
            loss_val = loss / entry_num
            return(loss_val.item())

    def train_total(self, epoch_num, patience):
        earlystopping = EarlyStopping(patience=patience, path=self.checkpoint)
        val_loss_list = deque(maxlen=patience) # val_loss_listをdequeとして初期化し、最大長さをpatienceに設定
        for epoch in range(epoch_num):
            loss = self.train_epoch()
            val_loss = self.evaluate(mode='validation')
            val_loss_list.append(val_loss)
            val_loss_mean = mean(val_loss_list)
            if self.use_val_loss_mean == True:
                earlystopping(val_loss_mean, self.scsurv)
            else:
                earlystopping(val_loss, self.scsurv)
            if earlystopping.early_stop:
                print(f"Early Stopping! at {epoch} epoch, best score={earlystopping.best_score}")
                break
            if self.scsurv.mode == 'beta_z':
                if epoch % 10 == 0:
                    total_data_loss = self.evaluate_train()
                    train_c_index = total_data_loss
                    val_c_index = val_loss
                    print(f'epoch {epoch}: train c-index {-train_c_index} validation c-index {-val_c_index}')
            elif epoch % 50 == 0:
                print(f'epoch {epoch}: train loss {loss} validation loss {val_loss}')
            if math.isnan(loss):
                print('loss is nan')
                break
        return epoch

    def initialize_optimizer(self, lr):
        self.scsurv_optimizer = torch.optim.AdamW(self.scsurv.parameters(), lr=lr)

    def initialize_loader(self, x_batch_size):
        self.x_data_manager.initialize_loader(x_batch_size)

    def calc_hazard_loss_beta_z(self, linear_predictor, survival_time, event_observed, bulk_idx, method):
        device = survival_time.device
        survival_time = survival_time[bulk_idx]
        linear_predictor = linear_predictor[bulk_idx].to(device)
        event_observed = event_observed[bulk_idx].to(device)
        hazard_ratio = torch.exp(linear_predictor)
        
        sorted_idx = torch.argsort(survival_time)
        sorted_survival_time = survival_time[sorted_idx]
        sorted_linear_predictor = linear_predictor[sorted_idx]
        sorted_hazard_ratio = hazard_ratio[sorted_idx]
        sorted_event_observed = event_observed[sorted_idx]
        
        log_risk = torch.zeros_like(sorted_survival_time)

        unique_times, counts = torch.unique_consecutive(sorted_survival_time[sorted_event_observed == 1], return_counts=True)
        if method == 'breslow':
            for t, n_events in zip(unique_times, counts):
                risk_set = sorted_survival_time >= t
                risk_sum = sorted_hazard_ratio[risk_set].sum()
                log_risk[(sorted_survival_time == t) & (sorted_event_observed == 1)] = torch.log(torch.clamp(risk_sum, min=1e-8))
        elif method == 'efron':
            for t, n_events in zip(unique_times, counts):
                risk_set = sorted_survival_time >= t
                tie_hazards = sorted_hazard_ratio[(sorted_survival_time == t) & (sorted_event_observed == 1)]
                risk_sum = sorted_hazard_ratio[risk_set].sum()
                tie_sum = tie_hazards.sum()
                if n_events == 1:
                    log_risk[(sorted_survival_time == t) & (sorted_event_observed == 1)] = torch.log(torch.clamp(risk_sum, min=1e-8))
                else:
                    efron_contributions = torch.stack([risk_sum - (j / n_events) * tie_sum 
                                                    for j in range(n_events)])
                    log_risk[(sorted_survival_time == t) & (sorted_event_observed == 1)] = torch.log(torch.clamp(efron_contributions, min=1e-8))
        else:
            raise ValueError("Method must be either 'breslow' or 'efron'")
        
        uncensored_likelihood = sorted_linear_predictor - log_risk
        censored_likelihood = uncensored_likelihood * sorted_event_observed
        neg_likelihood = -torch.sum(censored_likelihood)
        
        return neg_likelihood, sorted_survival_time, sorted_linear_predictor, sorted_event_observed
        
    def calc_kld(self, qz):
        kld = -0.5 * (1 + qz.scale.pow(2).log() - qz.loc.pow(2) - qz.scale.pow(2))
        # pz = Normal(torch.zeros(qz.loc.shape).to(self.device), torch.ones(qz.scale.shape).to(self.device))
        # kld = kl_divergence(qz, pz) #これと同じ
        return kld
    
    def calc_nb_loss(self, ld, norm_mat, theta, obs):
        ld = norm_mat * ld
        ld = ld + 1.0e-10
        theta = theta + 1.0e-10
        lp =  ld.log() - (theta).log()
        p_z = dist.NegativeBinomial(theta, logits=lp)
        l = - p_z.log_prob(obs)
        return l
    
    def calc_poisson_loss(self, ld, norm_mat, obs):
        p_z = dist.Poisson(ld * norm_mat + 1.0e-10)
        l = - p_z.log_prob(obs)
        return l