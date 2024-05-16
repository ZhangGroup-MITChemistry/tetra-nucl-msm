import numpy as np
import pandas as pd
try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit
import mdtraj
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 12})

"""Compute the relative diffusion between nucleosomes."""

os.makedirs('pictures', exist_ok=True)

nopbc_data = np.load('167bps_nucleosome_sea_single_nucleosome_coordinates_withoutpbc_xyz.npy', allow_pickle=True).item()
#print(type(nopbc_data)) # nopbc_data is a dictionary

timestep = 10 * unit.femtosecond # simulation timestep
delta_n_steps = 500 # the number of steps between two frames
delta_time_ns = (delta_n_steps * timestep).value_in_unit(unit.nanosecond) # the time between two frames in nanoseconds

box_length = 55 * unit.nanometer # cubic box length in unit nm
V = box_length**3 # cubic box volume in unit nm^3

# compute distance matrices
# for simplicity, as single nucleosomes are all identical, we only keep two single nucleosomes to save time and memory
distances = {}
data_keys = list(nopbc_data.keys()) # list including integers from 0 to 529
n_frames_list = []
for k in data_keys:
    traj = nopbc_data[k] # ndarray of shape (n_frames, n_nucls, 3)
    traj = traj[:, :8, :] # only keep 8 nucleosomes, including the tetra-nucleosome and 4 single nucleosomes
    # construct distance matrices
    n_frames = traj.shape[0]
    n_atoms = traj.shape[1]
    distances[k] = np.zeros((n_frames, n_atoms, n_atoms))
    rows, cols = np.triu_indices(n_atoms, k=1)
    r = traj[:, rows, :] - traj[:, cols, :]
    distances[k][:, rows, cols] = np.linalg.norm(r, axis=-1)
    distances[k][:, cols, rows] = distances[k][:, rows, cols] # symmetric
    n_frames_list.append(n_frames)
min_n_frames = min(n_frames_list)
t = np.arange(min_n_frames) * delta_time_ns

# track relative diffusion
pair_dict = {0: [[0, 4], [0, 5], [0, 6]], 
             1: [[1, 4], [1, 5], [1, 6]], 
             2: [[4, 5], [4, 6], [5, 6]]}
for i in range(3):
    # i is pair type index
    for p in pair_dict[i]:
        a1, a2 = p
        MSD = np.zeros(min_n_frames)
        for k in data_keys:
            r = distances[k][:min_n_frames, a1, a2] - distances[k][0, a1, a2]
            MSD += r**2
        MSD /= len(data_keys)
        plt.plot(t, MSD, label=f'nucleosomes {a1 + 1}-{a2 + 1}')
    plt.xlabel('t (ns)')
    plt.ylabel('MSD (nm^2)')
    plt.legend()
    if i == 0:
        plt.title('MSD between nucleosome 1 and single nucleosomes')
        plt.tight_layout()
        plt.savefig('pictures/MSD_nucl_1_single_nucl.pdf')
    if i == 1:
        plt.title('MSD between nucleosome 2 and single nucleosomes')
        plt.tight_layout()
        plt.savefig('pictures/MSD_nucl_2_single_nucl.pdf')
    if i == 2:
        plt.title('MSD between single nucleosomes')
        plt.tight_layout()
        plt.savefig('pictures/MSD_single_nucls.pdf')
    plt.close()



