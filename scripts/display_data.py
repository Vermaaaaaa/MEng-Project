import h5py
import os
import sys


path = os.path.dirname(os.path.abspath(__file__))

def disp(filepath, num: int):
    with h5py.File(filepath, 'r') as file:
        shape = file['imu_gyr'].shape
        print(shape)
        imu_gyr = file['imu_gyr'][:num]
        opt_quat = file['opt_quat'][:num]
        sampling_rate = file.attrs["sampling_rate"]
        dt = 1.0 / sampling_rate
        
        sample = {'imu_gyr': imu_gyr, 'opt_quat': opt_quat, 'dt': dt, 'sampling_rate': sampling_rate}

        print(sample)


        
if __name__ == "__main__":
    filename = sys.argv[1]
    
    try:
        num = int(sys.argv[2])
    except Exception as e:
        print(f"Usage: python3 display_data.py <file.hdf5> <num_rows_int>, Exception: {e}")
    
    filepath = os.path.join(path, filename)

    try:
        out = disp(filepath, num)
    except Exception as e:
        print(f"Failed on {filepath}: {e}")
