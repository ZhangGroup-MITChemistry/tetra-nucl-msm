import numpy as np
import pandas as pd
try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 12})

os.makedirs('pictures', exist_ok=True)

nopbc_data = np.load('167bps_nucleosome_sea_single_nucleosome_coordinates_withoutpbc_xyz.npy', allow_pickle=True).item()
#print(type(nopbc_data)) # nopbc_data is a dictionary

timestep = 10 * unit.femtosecond
delta_n_steps = 500 # the number of steps between two frames
delta_time = (delta_n_steps * timestep).value_in_unit(unit.nanosecond) # the time between two frames in nanoseconds

box_length = 55 # cubic box length in unit nm
V = box_length**3 # cubic box volume in unit nm^3

# compute distance matrices
distances = {}
data_keys = list(nopbc_data.keys())
for k in data_keys:
    traj = nopbc_data[k]
    #print(type(traj))
    #print(traj.shape)
    # construct distance matrices
    n_frames = traj.shape[0]
    n_atoms = traj.shape[1]
    distances[k] = np.zeros((n_frames, n_atoms, n_atoms))
    rows, cols = np.triu_indices(n_atoms, k=1)
    r = traj[:, rows, :] - traj[:, cols, :]
    r = (r / box_length) % 1
    r[r >= 0.5] -= 1
    assert np.amin(r) >= -0.5
    assert np.amax(r) < 0.5
    r *= box_length
    distances[k][:, rows, cols] = np.linalg.norm(r, axis=-1)
    distances[k][:, cols, rows] = distances[k][:, rows, cols] # symmetric

# compute and draw radial distribution functions for each stage
stages = {0: [0, 20], 1: [20, 40], 2: [40, 50], 3: [60, 70]}
nucl_types_to_indices = {0: [0, 3], 1: [1, 2], 2: list(range(4, 30))}
hists = [[], [], [], []] # hists[i][j] is the histogram of stage i and reference nucleosome type j
assert len(hists) == len(stages)
bins = np.linspace(0.1, 20, 50) # start from 0.1 to avoid including reference nucleosome itself
x = (bins[1:] + bins[:-1]) / 2
for i in range(len(stages)):
    # i is the stage index
    start_frame_index = int(stages[i][0] / delta_time)
    end_frame_index = int(stages[i][1] / delta_time)
    print(f'Stage {i} uses frame {start_frame_index}-{end_frame_index}')
    for j in range(3):
        # j is the reference nucleosome type index
        samples = [] # collect all the samples with given reference nucleosome type and stage
        ref_nucl_indices = nucl_types_to_indices[j]
        for k in data_keys:
            traj = distances[k]
            n_frames = traj.shape[0]
            end_frame_index = min(n_frames - 1, end_frame_index)
            assert end_frame_index > start_frame_index
            mask = np.full((n_atoms, n_atoms), False)
            mask[ref_nucl_indices, :] = True
            np.fill_diagonal(mask, False)
            samples += traj[start_frame_index:(end_frame_index + 1), mask].flatten().tolist()
        samples = np.array(samples)
        hists[i].append(np.histogram(samples, bins=bins, density=True)[0])

rho = 30 / V
for j in range(3):
    # j is the reference nucleosome type index
    for i in range(len(stages)):
        # i is the stage index
        g = (30 - 1) * hists[i][j] / (4 * np.pi * rho * x**2)
        plt.plot(x, g, label=f'{stages[i][0]}-{stages[i][1]} ns')
    plt.legend()
    plt.xlabel('r (nm)')
    plt.ylabel('g(r)')
    if j == 0:
        plt.title('Nucleosome 1, 4 as references')
    elif j == 1:
        plt.title('Nucleosome 2, 3 as references')
    else:
        plt.title('Single nucleosomes as references')
    plt.savefig(f'pictures/radial_distribution_type_{j}.pdf')
    plt.close()





