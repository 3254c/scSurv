import torch
import torch.distributions as dist
import torch.nn as nn
from torch.nn.parameter import Parameter
from torch.nn import init

class LinearGELU(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(LinearGELU, self).__init__()
        self.f = nn.Sequential(
            nn.Linear(input_dim, output_dim),
            nn.LayerNorm(output_dim, elementwise_affine=False),
            nn.GELU())

    def forward(self, x):
        h = self.f(x)
        return(h)

class SeqNN(nn.Module):
    def __init__(self, num_steps, dim):
        super(SeqNN, self).__init__()
        modules = [
            LinearGELU(dim, dim)
            for _ in range(num_steps)
        ]
        self.f = nn.Sequential(*modules)

    def forward(self, pre_h):
        post_h = self.f(pre_h)
        return(post_h)


class Encoder(nn.Module):
    def __init__(self, num_h_layers, x_dim, h_dim, z_dim, enc_dist=dist.Normal):
        super(Encoder, self).__init__()
        self.x2h = LinearGELU(x_dim, h_dim)
        self.seq_nn = SeqNN(num_h_layers - 1, h_dim)
        self.h2mu = nn.Linear(h_dim, z_dim)
        self.h2logvar = nn.Linear(h_dim, z_dim)
        self.softplus = nn.Softplus()
        self.dist = enc_dist

    def forward(self, x):
        pre_h = self.x2h(x)
        post_h = self.seq_nn(pre_h)
        mu = self.h2mu(post_h)
        logvar = self.h2logvar(post_h) #対数分散 log(sigma^2)
        sigma = self.softplus(logvar)
        qz = self.dist(mu, sigma)
        z = qz.rsample()
        return(z, qz)

class Decoder_z(nn.Module):
    def __init__(self, num_h_layers, z_dim, h_dim, x_dim):
        super(Decoder_z, self).__init__()
        self.z2h = LinearGELU(z_dim, h_dim)
        self.seq_nn = SeqNN(num_h_layers - 1, h_dim)
        self.h2ld = nn.Linear(h_dim, x_dim)
        self.softplus = nn.Softplus()

    def forward(self, z):
        pre_h = self.z2h(z)
        post_h = self.seq_nn(pre_h)
        ld = self.h2ld(post_h)
        softplus_ld = self.softplus(ld)
        return(softplus_ld)

class Decoder_p_softmax(nn.Module):
    def __init__(self, num_h_layers, z_dim, h_dim, x_dim):
        super(Decoder_p_softmax, self).__init__()
        self.z2h = LinearGELU(z_dim, h_dim)
        self.seq_nn = SeqNN(num_h_layers - 1, h_dim)
        self.h2ld = nn.Linear(h_dim, x_dim)
        self.softmax = nn.Softmax(dim=1)
        self.softplus = nn.Softplus()
    
    def forward(self, z):
        pre_h = self.z2h(z)
        post_h = self.seq_nn(pre_h)
        ld = self.h2ld(post_h)
        ld_t = ld.transpose(0, 1)
        p_softmax = self.softmax(ld_t)
        return(p_softmax)
    
class Decoder_beta_z(nn.Module):
    def __init__(self, num_h_layers, z_dim, h_dim, x_dim):
        super(Decoder_beta_z, self).__init__()
        self.z2h = LinearGELU(z_dim, h_dim)
        self.seq_nn = SeqNN(num_h_layers - 1, h_dim)
        self.h2ld = nn.Linear(h_dim, x_dim)
        self.softplus = nn.Softplus()

    def forward(self, z):
        pre_h = self.z2h(z)
        post_h = self.seq_nn(pre_h)
        ld = self.h2ld(post_h)
        return(ld)

class scSurv(nn.Module):
    def __init__(self, bulk_num, spatial_num, batch_onehot_dim, x_dim, z_dim, h_dim, num_enc_z_layers, num_dec_z_layers, num_dec_p_layers, num_dec_b_layers):
        super(scSurv, self).__init__()
        self.bulk_num = bulk_num
        self.spatial_num = spatial_num
        self.enc_z = Encoder(num_enc_z_layers, x_dim + batch_onehot_dim, h_dim, z_dim)
        self.dec_z2x = Decoder_z(num_dec_z_layers, z_dim + batch_onehot_dim, h_dim, x_dim)
        self.dec_z2p_bulk = Decoder_p_softmax(num_dec_p_layers, z_dim, h_dim, bulk_num)
        self.dec_z2p_spatial = Decoder_p_softmax(num_dec_p_layers, z_dim, h_dim, spatial_num)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.softplus = nn.Softplus()
        self.num_dec_b_layers = num_dec_b_layers
        self.z_dim = z_dim
        self.h_dim = h_dim
        self.log_spatial_coeff = Parameter(torch.Tensor(x_dim))
        self.log_spatial_coeff_add = Parameter(torch.Tensor(x_dim))
        self.log_spatial_theta=  Parameter(torch.Tensor(x_dim))
        self.log_bulk_coeff = Parameter(torch.Tensor(x_dim))
        self.log_bulk_coeff_add = Parameter(torch.Tensor(x_dim))
        self.log_bulk_theta=  Parameter(torch.Tensor(x_dim))
        self.logtheta_x =  Parameter(torch.Tensor(x_dim))
        self.mode = 'sc'
        self.dec_beta_z = Decoder_beta_z(num_dec_b_layers, z_dim, h_dim, 1)
        self.reset_parameters_zeros()

    def reset_parameters_zeros(self):
        init.constant_(self.log_spatial_coeff, 0)
        init.constant_(self.log_spatial_coeff_add, 0)
        init.constant_(self.log_spatial_theta, 0)
        init.constant_(self.log_bulk_coeff, 0)
        init.constant_(self.log_bulk_coeff_add, 0)
        init.constant_(self.log_bulk_theta, 0)
        init.constant_(self.logtheta_x, 0)

    def forward(self, x, batch_onehot):
        xb = torch.cat([x, batch_onehot], dim=-1)
        z, qz = self.enc_z(xb)
        zb = torch.cat([z, batch_onehot], dim=-1)
        x_hat = self.dec_z2x(zb)
        p_bulk = self.dec_z2p_bulk(z)
        bulk_coeff = self.softplus(self.log_bulk_coeff) #gene size
        bulk_coeff_add = self.softplus(self.log_bulk_coeff_add) #bulk
        bulk_hat = torch.matmul(p_bulk, x_hat * bulk_coeff) + bulk_coeff_add #spatialの発現のmean parameter #matmulは行列の積
        if self.spatial_num == 0:
            p_spatial = None
            spatial_coeff = None
            spatial_coeff_add = None
            spatial_hat = None
        else:
            p_spatial = self.dec_z2p_spatial(z) # * (z.shape[0] * self.x_num) #scRNAseq中の各細胞がbulkに占める割合 bulkごと和が1になる #sc_adata.obsm['p_mat'].sum(axis=0)は(1, 1, ..., 1)bulkデータ次元 つまりbulkデータ内の各細胞の割合
            spatial_coeff = self.softplus(self.log_spatial_coeff) #gene size
            spatial_coeff_add = self.softplus(self.log_spatial_coeff_add) #spatial
            spatial_hat = torch.matmul(p_spatial, x_hat * spatial_coeff) + spatial_coeff_add #spatialの発現のmean parameter #matmulは行列の積
        theta_x = self.softplus(self.logtheta_x)
        theta_bulk = self.softplus(self.log_bulk_theta)
        theta_spatial = self.softplus(self.log_spatial_theta)
        beta_z = self.dec_beta_z(z).squeeze()
        return(z, qz, x_hat, p_bulk, p_spatial, bulk_hat, spatial_hat, theta_x, theta_bulk, theta_spatial, beta_z)

    def sc_mode(self):
        self.mode = 'sc'
        for parameter in self.parameters():
            parameter.requires_grad = False
        for parameter in self.enc_z.parameters():
            parameter.requires_grad = True
        for parameter in self.dec_z2x.parameters():
            parameter.requires_grad = True
        self.logtheta_x.requires_grad = True

    def bulk_mode(self):
        self.mode = 'bulk'
        for parameter in self.parameters():
            parameter.requires_grad = False
        for parameter in self.dec_z2p_bulk.parameters():
            parameter.requires_grad = True
        self.log_bulk_coeff.requires_grad = True
        self.log_bulk_coeff_add.requires_grad = True
        self.log_bulk_theta.requires_grad = True

    def spatial_mode(self):
        self.mode = 'spatial'
        for parameter in self.parameters():
            parameter.requires_grad = False
        for parameter in self.dec_z2p_spatial.parameters():
            parameter.requires_grad = True
        self.log_spatial_coeff.requires_grad = True
        self.log_spatial_coeff_add.requires_grad = True
        self.log_spatial_theta.requires_grad = True
        
    def hazard_beta_z_mode(self):
        self.mode = 'beta_z'
        for parameter in self.parameters():
            parameter.requires_grad = False
        for parameter in self.dec_beta_z.parameters():
            parameter.requires_grad = True