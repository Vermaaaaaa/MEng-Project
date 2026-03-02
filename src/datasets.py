import torch
from torch.utils.data import Dataset
import h5py


# ''' TODO ''' 
# Normalise the data

class BroadDataset(Dataset):
    def __init__(self, filepaths, mode, window, stride=None, dtype=torch.float32, transform=None):
        super().__init__()
        
        self.mode = mode
        self.dtype = dtype
        self.transform = transform
        self.window = int(window)
        self.filepaths = list(filepaths)
        if stride is None:
            self.stride = int(window)
        else:
            self.stride = int(stride)

        self._handles = {}
        self.index = []

        self.noise_x = 0.10
        self.noise_y = 0.09
        self.noise_z = 0.12


        with h5py.File(filepaths[0], 'r') as file:
            self.sampling_rate = file.attrs['sampling_rate']
        
        self.dt = 1.0 / self.sampling_rate


        for filepath in filepaths:
            with h5py.File(filepath, 'r') as f:
                length = f['imu_gyr'].shape[0]
                num_windows = length // window
                for i in range(0,num_windows):
                    self.index.append((filepath, i, length))
        

        

    def __len__(self):
        return len(self.index)
    
    def _get_file_handle(self, filepath):
        if filepath not in self._handles:
            self._handles[filepath] = h5py.File(filepath, 'r')

        return self._handles[filepath]
        
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        if self.mode == 'train':
            filepath, window_index, length = self.index[idx]
            file = self._get_file_handle(filepath)
            
            max_start = length - self.window
            start = torch.randint(0, max_start + 1, (1, )).item()
            end = start + self.window
            imu_gyr = file['imu_gyr'][start:end]
            opt_quat = file['opt_quat'][start:end]

            noise = self.generate_noise()

            imu_gyr = torch.as_tensor(imu_gyr, dtype=self.dtype)
            opt_quat = torch.as_tensor(opt_quat, dtype=self.dtype)
            imu_gyr = imu_gyr + noise

        else:
            filepath, window_index, length = self.index[idx]
            file = self._get_file_handle(filepath)

            start = window_index * self.window

            imu_gyr = file['imu_gyr'][start: start+self.window]
            opt_quat = file['opt_quat'][start: start+self.window]

            imu_gyr = torch.as_tensor(imu_gyr, dtype=self.dtype)
            opt_quat = torch.as_tensor(opt_quat, dtype=self.dtype)


        sample = {
            'imu_gyr': imu_gyr, 
            'opt_quat': opt_quat, 
            'dt': self.dt, 
            'sampling_rate': self.sampling_rate
            }

        if self.transform:
           sample = self.transform(sample)
        return sample
    
    def generate_noise(self):
        std = torch.tensor([self.noise_x, self.noise_y, self.noise_z]).expand(self.window, 3)
        noise = torch.normal(mean=0., std=std)

        return noise

    def __del__(self):
        for file in self._handles.values():
            try:
                file.close()
            except Exception:
                pass
        self._handles.clear()
    
    
