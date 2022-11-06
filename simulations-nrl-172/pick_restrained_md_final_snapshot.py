import numpy as np
import mdtraj as md
import os

if not os.path.exists('restrained-md-final-snapshots'):
    os.makedirs('restrained-md-final-snapshots')

# 4643 simulations in all
mdconvert_path = '/home/xclin/bin/anaconda2/bin/mdconvert'
for i in range(4643):
    output_dcd_path = f'restrained-md-final-snapshots/snapshot{i}.dcd'
    if os.path.exists(output_dcd_path):
        continue
    dcd_path = 'restrained-md/run%05d/sim/run000/DUMP_FILE.dcd' % i
    if not os.path.exists(dcd_path):
        print(f'{dcd_path} does not exist!')
        continue
    traj = md.load_dcd(dcd_path, 'chromatin-172x4.pdb')
    n_frames = traj.n_frames
    if n_frames < 11:
        print(f'Warning: {dcd_path} does not have enough frames!')
        print(f'Do not pick the final snapshot from {dcd_path}')
        continue
    cmd = f'{mdconvert_path} -o {output_dcd_path} -i -1 {dcd_path}'
    os.system(cmd)
    
        



