import numpy as np
import pandas as pd
import mdtraj

"""
Prepare the pdb file for drawing the condensate system
"""

traj = mdtraj.load_pdb('condensate-traj-run382/start.pdb')
n_frames = traj.n_frames
a = 55.0 # cubic box length in unit nm
traj.unitcell_lengths = np.array([[a, a, a]] * n_frames)
traj.unitcell_angles = np.array([[90.0, 90.0, 90.0]] * n_frames)
traj.save_pdb('condensate-traj-run382/start_with_box.pdb')


