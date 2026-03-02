import h5py
import numpy as np
import os

required = {'imu_gyr', 'imu_acc', 'imu_mag', 'opt_quat', 'movement'}

path = os.path.dirname(os.path.abspath(__file__))

def clean_broad_data(filepath):
    with h5py.File(filepath, 'r') as f: 
        missing = required - set(f.keys())
        if missing:
            raise KeyError(f"Missing datasets: {missing}")


        opt_quat = f['opt_quat']
        movement = f['movement']
        imu_gyr = f['imu_gyr']
        imu_acc = f['imu_acc']
        imu_mag = f['imu_mag']
        sampling_rate = f.attrs['sampling_rate']
        

        mask = ~np.isnan(opt_quat).any(axis=1)

        cleaned_opt_quat = opt_quat[mask]
        cleaned_movement = movement[mask]
        cleaned_imu_gyr = imu_gyr[mask]
        cleaned_imu_acc = imu_acc[mask]
        cleaned_imu_mag = imu_mag[mask]

        print("Removed rows:", opt_quat.shape[0] - cleaned_opt_quat.shape[0])


        if (cleaned_opt_quat.shape[0] 
            == cleaned_imu_gyr.shape[0] 
            == cleaned_imu_acc.shape[0] 
            == cleaned_imu_mag.shape[0] 
            == cleaned_movement.shape[0] 
            and not np.isnan(cleaned_opt_quat).any()):

            new_dir = os.path.join(path, 'broad_cleaned')
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            
            base = os.path.basename(filepath)
            name, ext = os.path.splitext(base)
            new_filepath = os.path.join(new_dir, f"{name}_cleaned{ext}")

            with h5py.File(new_filepath, 'w') as file:
                file.create_dataset('imu_gyr', data=cleaned_imu_gyr, compression='gzip', compression_opts=9)
                file.create_dataset('imu_acc', data=cleaned_imu_acc, compression='gzip', compression_opts=9)
                file.create_dataset('imu_mag', data=cleaned_imu_mag, compression='gzip', compression_opts=9)
                file.create_dataset('opt_quat', data=cleaned_opt_quat, compression='gzip', compression_opts=9)
                file.create_dataset('movement', data=cleaned_movement, compression='gzip', compression_opts=9)
                file.attrs['sampling_rate'] = sampling_rate
                file.attrs['reference_frame'] = 'ENU'

            with h5py.File(filepath, 'r') as f1, h5py.File(new_filepath, 'r') as f2:
                opt_quat_original = f1['opt_quat']
                mask = ~np.isnan(opt_quat_original).any(axis=1)
                cleaned_opt_quat_original = opt_quat_original[mask]

                cleaned_opt_quat_test = f2['opt_quat']

                if not np.allclose(cleaned_opt_quat_test[:], cleaned_opt_quat_original[:]):
                        raise AssertionError("Values are not the same")

    return new_filepath




if __name__ == "__main__":
    for filename in os.listdir(path):
        if not filename.endswith(".hdf5"):
            continue

        if filename.endswith("_cleaned.hdf5"):
            continue

        filepath = os.path.join(path, filename)

        try:
            print(f"Cleaning {filename}...")
            out = clean_broad_data(filepath)
            print(f"Saved → {out}")
        except Exception as e:
            print(f"Failed on {filename}: {e}")



