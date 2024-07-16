import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 12})

"""
Compute the PMF along specific CVs by reweighting the trajectories. 
The trajectories are many short unbiased simulations. 
Following Yunrui's code.
"""

output_dir = 'reweight-results'
os.makedirs(output_dir, exist_ok=True)

## load the six inter-nucleosome distances based trajectories
# note the 6 distances are arranged as d12, d13, d14, d23, d24, d34
distrajs = np.load("./177bps_6pdis_mirror_cut_trajs.npy", allow_pickle=True).item()
distrajs = [value for value in distrajs.values()]
distrajs = np.concatenate(distrajs)
#print(distrajs.shape)

## load the equilibrium population of microstate-MSM and the microstate-assigned-trajectories
stationary_population = np.loadtxt("./177bps_0.25ns3tics_1000kmeans_equilibrium_population.txt")
stationary_population /= np.sum(stationary_population) # ensure normalized
ctrajs = np.load("./177bps_0.25ns3tics_1000kmeans_trajs.npy", allow_pickle=True)
ctrajs = np.concatenate(ctrajs).astype(int) # ctrajs is the microstate index of each snapshot
num_states = np.max(ctrajs) + 1
#print(f'{num_states} microstates in all')

## randomly regenerate index trajectories to match the size of the microstate-assigned trajectories, 
## ensuring the frequency of each microstate is proportional to its corresponding equilibrium population.
## the number of microstates for NRL=177 system is 1000
#n_frames = int(10**6) # the number of frames used for computing PMF
n_frames = ctrajs.shape[0] # the target total number of frames for analysis
selected_distrajs = []
for i in range(num_states):
    n_i = int(round(n_frames * stationary_population[i]))
    distrajs_i = distrajs[ctrajs == i]
    selection = np.random.choice(len(distrajs_i), n_i)
    selected_distrajs.append(distrajs_i[selection])
selected_distrajs = np.concatenate(selected_distrajs)

# compute PMF
# for d13
d13_samples = selected_distrajs[:, 1]
p, bin_edges = np.histogram(d13_samples, bins=100, density=True)
r = (bin_edges[:-1] + bin_edges[1:]) / 2
pmf = -np.log(p)
pmf -= np.min(pmf) # set the minimum to 0
df_d13_pmf = pd.DataFrame({'r (nm)': r, 'PMF (kT)': pmf})
df_d13_pmf.to_csv(f'{output_dir}/d13_pmf.csv', index=False)

with PdfPages(f'{output_dir}/d13.pdf') as pdf:
    plt.plot(r, p)
    plt.xlabel('d13 (nm)')
    plt.ylabel('Probability density')
    plt.title(f'Probability density along d13')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.plot(r, pmf)
    plt.xlabel('d13 (nm)')
    plt.ylabel('PMF (kT)')
    plt.title(f'PMF along d13')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

# for d23
d23_samples = selected_distrajs[:, 3]
p, bin_edges = np.histogram(d23_samples, bins=100, density=True)
r = (bin_edges[:-1] + bin_edges[1:]) / 2
pmf = -np.log(p)
pmf -= np.min(pmf) # set the minimum to 0
df_d23_pmf = pd.DataFrame({'r (nm)': r, 'PMF (kT)': pmf})
df_d23_pmf.to_csv(f'{output_dir}/d23_pmf.csv', index=False)

with PdfPages(f'{output_dir}/d23.pdf') as pdf:
    plt.plot(r, p)
    plt.xlabel('d23 (nm)')
    plt.ylabel('Probability density')
    plt.title(f'Probability density along d23')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.plot(r, pmf)
    plt.xlabel('d23 (nm)')
    plt.ylabel('PMF (kT)')
    plt.title(f'PMF along d23')
    plt.tight_layout()
    pdf.savefig()
    plt.close()



