import numpy as np
import os

dist_data = []
for i in range(10000):
    colvar_path = 'plumed-driver-output/COLVAR%05d' % i
    if os.path.exists(colvar_path):
        colvar = np.loadtxt(colvar_path, skiprows=1)
        if colvar.ndim == 1:
            d13, d24 = colvar[2], colvar[5]
            if d13 >= d24:
                dist_data.append(colvar[1:])
        else:
            print(f'check the shape of {colvar_path}')
    else:
        print(f'{colvar_path} does not exist')

dist_data = np.array(dist_data)
n_selected_samples = dist_data.shape[0]
print(f'{n_selected_samples} selected samples with have d13 >= d24')
print('lengths in unit nm')
header = 'd12 d13 d14 d23 d24 d34'
np.savetxt('selected_sample_dist.txt', dist_data, fmt='%.6f', header=header)

