import torch
import numpy as np

class ScDataSet(torch.utils.data.Dataset):
    def __init__(self, x, xnorm_mat, batch_onehot):
        self.x = x
        self.xnorm_mat = xnorm_mat
        self.batch_onehot = batch_onehot

    def __len__(self):
        return(self.x.shape[0])

    def __getitem__(self, idx):
        idx_x = self.x[idx]
        idx_xnorm_mat = self.xnorm_mat[idx]
        idx_batch_onehot = self.batch_onehot[idx]
        return(idx_x, idx_xnorm_mat, idx_batch_onehot)

class ScDataManager():
    def __init__(self, x_count, batch_size, batch_onehot):
        validation_ratio = 0.1
        test_ratio = 0.05

        x_count = x_count.float()
        self.batch_onehot = batch_onehot
        xnorm_mat = torch.mean(x_count, dim=1).view(-1, 1)
        total_num = x_count.size()[0]
        self.x_count = x_count
        self.idx_init = torch.tensor(np.array(range(total_num)))
        self.xnorm_mat = xnorm_mat

        validation_num = int(total_num * validation_ratio)
        test_num = int(total_num * test_ratio)

        np.random.seed(42)
        idx = np.random.permutation(np.arange(total_num))
        self.idx = torch.tensor(idx)

        validation_idx, test_idx, train_idx = idx[:validation_num], idx[validation_num:(validation_num +  test_num)], idx[(validation_num +  test_num):]
        self.validation_idx, self.test_idx, self.train_idx = torch.tensor(validation_idx), torch.tensor(test_idx), torch.tensor(train_idx)

        self.train_x = x_count[train_idx]
        self.train_xnorm_mat = xnorm_mat[train_idx]
        self.validation_x = x_count[validation_idx]
        self.validation_xnorm_mat = xnorm_mat[validation_idx]
        self.test_x = x_count[test_idx]
        self.test_xnorm_mat = xnorm_mat[test_idx]
        self.train_batch_onehot = batch_onehot[train_idx]
        self.validation_batch_onehot = batch_onehot[validation_idx]
        self.test_batch_onehot = batch_onehot[test_idx]
        self.train_eds = ScDataSet(x_count[train_idx], xnorm_mat[train_idx], batch_onehot[train_idx])
        self.train_loader = torch.utils.data.DataLoader(self.train_eds, batch_size=batch_size, shuffle=True, num_workers=8, drop_last=True, pin_memory=True)

    def initialize_loader(self, batch_size):
        self.train_loader = torch.utils.data.DataLoader(self.train_eds, batch_size=batch_size, shuffle=True, num_workers=8, drop_last=True, pin_memory=True)

class BulkDataManager():
    def __init__(self, bulk_count, survival_time, cutting_off_0_1):
        bulk_count = bulk_count.float()
        survival_time = survival_time.float()
        bnorm_mat = torch.mean(bulk_count, dim=1).view(-1, 1)
        self.bulk_count = bulk_count
        self.survival_time = survival_time
        self.cutting_off_0_1 = cutting_off_0_1
        self.bulk_norm_mat = bnorm_mat
        self.train_idx = torch.tensor([])
        self.validation_idx = torch.tensor([])
        self.test_idx = torch.tensor([])
    
    def bulk_split(self, bulk_seed, bulk_validation_num_or_ratio, bulk_test_num_or_ratio, censor_np):
        idx_0 = np.where(censor_np == 0)[0]
        idx_1 = np.where(censor_np == 1)[0]
        N0 = len(idx_0)
        N1 = len(idx_1)
        total_num = N0 + N1
        if bulk_validation_num_or_ratio < 1:
            validation_ratio = bulk_validation_num_or_ratio
        else:
            validation_ratio = bulk_validation_num_or_ratio / total_num
        if bulk_test_num_or_ratio < 1:
            test_ratio = bulk_test_num_or_ratio
        else:
            test_ratio = bulk_test_num_or_ratio / total_num
        
        np.random.seed(bulk_seed)

        idx_0_permuted = np.random.permutation(idx_0)
        N0_val = int(N0 * validation_ratio)
        N0_test = int(N0 * test_ratio)
        idx_0_validation = idx_0_permuted[:N0_val]
        idx_0_test = idx_0_permuted[N0_val:N0_val + N0_test]
        idx_0_train = idx_0_permuted[N0_val + N0_test:]

        idx_1_permuted = np.random.permutation(idx_1)
        N1_val = int(N1 * validation_ratio)
        N1_test = int(N1 * test_ratio)
        idx_1_validation = idx_1_permuted[:N1_val]
        idx_1_test = idx_1_permuted[N1_val:N1_val + N1_test]
        idx_1_train = idx_1_permuted[N1_val + N1_test:]

        validation_idx = np.concatenate([idx_0_validation, idx_1_validation])
        test_idx = np.concatenate([idx_0_test, idx_1_test])
        train_idx = np.concatenate([idx_0_train, idx_1_train])

        np.random.shuffle(validation_idx)
        np.random.shuffle(test_idx)
        np.random.shuffle(train_idx)

        self.train_idx = torch.tensor(train_idx)
        self.validation_idx = torch.tensor(validation_idx)
        self.test_idx = torch.tensor(test_idx)
            
class SpatialDataManager():
    def __init__(self, spatial_count):
        spatial_count = spatial_count.float()
        spatial_norm_mat = torch.mean(spatial_count, dim=1).view(-1, 1)
        self.spatial_count = spatial_count
        self.spatial_norm_mat = spatial_norm_mat