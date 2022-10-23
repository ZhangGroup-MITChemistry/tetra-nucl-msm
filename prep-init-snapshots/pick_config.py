import numpy as np
import sys
import os
import shutil

mdconvert_path = '/home/xclin/bin/anaconda2/bin/mdconvert'
plumed_path = '/home/xclin/bin/plumed246_serial/bin/plumed'

if not os.path.exists('all-snapshots'):
    os.makedirs('all-snapshots')
if not os.path.exists('selected-snapshots'):
    os.makedirs('selected-snapshots')
if not os.path.exists('plumed-driver-output'):
    os.makedirs('plumed-driver-output')

# pick all the configurations from previous tetranucleosome restraint MD simulations
# then run plumed driver to check the CV values
# finally select configurations with d13 >= d24
main_dir = '/nobackup1/binz/chromatin_fiber_CG/tetra/nnet/sim6_ppInt'
selected_snapshot_index = 0
with open('template_plumed.txt', 'r') as template_plumed:
    template_plumed_text = template_plumed.read()
for i in range(10000):
    input_dcd_path = '%s/run%05d/DUMP_FILE.dcd' % (main_dir, i)
    output_dcd_path = 'all-snapshots/snapshot_%d.dcd' % i
    if os.path.exists(input_dcd_path):
        cmd = '%s -o %s -f -i -1 %s' % (mdconvert_path, output_dcd_path, input_dcd_path)
        os.system(cmd)
        colvar_path = 'plumed-driver-output/COLVAR%05d' % i
        with open('plumed.txt', 'w') as plumed:
            plumed_text = template_plumed_text.replace('OUTPUT_COLVAR_PATH', colvar_path)
            plumed.write(plumed_text)
        cmd = '%s driver --plumed plumed.txt --mf_dcd %s' % (plumed_path, output_dcd_path)
        os.system(cmd)
        colvar = np.loadtxt(colvar_path, skiprows=1)
        d13 = colvar[2]
        d24 = colvar[5]
        if d13 >= d24:
            shutil.copyfile(output_dcd_path, 'selected-snapshots/snapshot_%d.dcd' % selected_snapshot_index)
            selected_snapshot_index += 1
    else:
        print('%s does not exist' % dcd_path)

