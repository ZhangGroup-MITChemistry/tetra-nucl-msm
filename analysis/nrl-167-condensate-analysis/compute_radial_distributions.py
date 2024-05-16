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
from openabc.lib import NA

"""
Compute the radial distributions of nucleosomes.
"""

os.makedirs('pictures', exist_ok=True)

nopbc_data = np.load('167bps_nucleosome_sea_single_nucleosome_coordinates_withoutpbc_xyz.npy', allow_pickle=True).item()
#print(type(nopbc_data)) # nopbc_data is a dictionary

timestep = 10 * unit.femtosecond # simulation timestep
delta_n_steps = 500 # the number of steps between two frames
delta_time_ns = (delta_n_steps * timestep).value_in_unit(unit.nanosecond) # the time between two frames in nanoseconds

box_length = 55 * unit.nanometer # cubic box length in unit nm
V = box_length**3 # cubic box volume in unit nm^3

# compute distance matrices
distances = {}
data_keys = list(nopbc_data.keys()) # list including integers from 0 to 529
for k in data_keys:
    traj = nopbc_data[k] # ndarray of shape (n_frames, n_nucls, 3), and n_nucls = 30
    # construct distance matrices
    n_frames = traj.shape[0]
    n_atoms = traj.shape[1]
    assert n_atoms == 30
    distances[k] = np.zeros((n_frames, n_atoms, n_atoms))
    rows, cols = np.triu_indices(n_atoms, k=1)
    r = traj[:, rows, :] - traj[:, cols, :]
    L = box_length.value_in_unit(unit.nanometer)
    r = (np.mod((r / L) + 0.5, 1.0) - 0.5) * L
    distances[k][:, rows, cols] = np.linalg.norm(r, axis=-1)
    distances[k][:, cols, rows] = distances[k][:, rows, cols] # symmetric

# compute and draw radial distribution functions for each stage
# divide into 4 stages based on the time intervals
# divide nucleosomes into 3 types based on topology
stages = {0: [0, 20], 1: [20, 40], 2: [40, 50], 3: [60, 70]} # in time unit ns
nucl_type_to_ref_index = {0: 0, 1: 1, 2: 4} # use nucleosome as the reference for each type
hists = {}
bins = np.linspace(0.1, 20, 50) # start from 0.1 to avoid including reference nucleosome itself
x = (bins[1:] + bins[:-1]) / 2
for i in range(len(stages)):
    # i is the stage index
    hists[i] = {}
    start_frame_index = int(stages[i][0] / delta_time_ns)
    end_frame_index = int(stages[i][1] / delta_time_ns)
    print(f'Stage {i} uses frame {start_frame_index}-{end_frame_index}')
    for j in range(3):
        # j is the reference nucleosome type index
        samples = [] # collect all the samples with given reference nucleosome type and stage
        ref_nucl_index = nucl_type_to_ref_index[j]
        for k in data_keys:
            d = distances[k]
            n_frames = d.shape[0]
            end_frame_index = min(n_frames - 1, end_frame_index)
            assert end_frame_index > start_frame_index
            mask = np.full((n_atoms, n_atoms), False)
            mask[ref_nucl_index, :] = True
            np.fill_diagonal(mask, False)
            samples += d[start_frame_index:(end_frame_index + 1), mask].flatten().tolist()
        samples = np.array(samples)
        hists[i][j] = np.histogram(samples, bins=bins, density=True)[0]

rho = (30 / (NA * V)).value_in_unit(unit.millimolar) # rho in unit mM
for j in range(3):
    # j is the reference nucleosome type index
    ref_nucl_index = nucl_type_to_ref_index[j]
    for i in range(4):
        # i is the stage index
        g = (30 - 1) * hists[i][j] / (4 * np.pi * rho * x**2)
        plt.plot(x, g, label=f'{stages[i][0]}-{stages[i][1]} ns')
    plt.legend()
    plt.xlabel('r (nm)')
    plt.ylabel('g(r) (mM)')
    plt.title(f'Nucleosome {ref_nucl_index + 1} as reference')
    plt.tight_layout()
    plt.savefig(f'pictures/radial_distribution_nucl_{ref_nucl_index + 1}_ref.pdf')
    plt.close()


