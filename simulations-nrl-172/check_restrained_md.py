import numpy as np
import sys
import os
import glob

# 4643 sets of simulations in all
for i in range(4643):
    dcd_path = 'restrained-md/run%05d/sim/run000' % i
    if not os.path.exists(dcd_path):
        print(f'{dcd_path} does not exist!')



