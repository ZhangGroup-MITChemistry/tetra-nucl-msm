import numpy as np
import sys
import os
import math
import mdtraj

n_jobs = 4643
n_data_dir = math.ceil(n_jobs/500)

data_dir_index = int(sys.argv[1])

flag = True
for i in range(500*data_dir_index, 500*(data_dir_index + 1)):
    if i < n_jobs:
        sim_dir = 'unbiased-md/unbiased-md-data-%d/run%05d' % (data_dir_index, i)
        if os.path.exists(sim_dir):
            for j in range(2):
                dcd = '%s/sim/run%03d/DUMP_FILE.dcd' % (sim_dir, j)
                if os.path.exists(dcd):
                    traj = mdtraj.load_dcd(dcd, top='chromatin-172x4.pdb')
                    if traj.n_frames == 501:
                        pass
                    else:
                        flag = False
                        print(f'Snapshot number in {dcd} is incorrect')
                else:
                    flag = False
                    print(f'{dcd} does not exist')
        else:
            flag = False
            print(f'{sim_dir} does not exist')

if flag:
    print(f'data in unbiased-md/unbiased-md-data-{data_dir_index} are ready for sharing')

