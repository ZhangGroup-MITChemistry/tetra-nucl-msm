import numpy as np
import pandas as pd
import shutil
import sys
import os

sys.path.append('/home/gridsan/sliu/Projects/smog-3spn2-openmm')
import OpenSMOG3SPN2.utils.helper_functions as helper_functions

amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
               'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
               'LEU', 'LYS', 'MET', 'PHE', 'PRO',
               'SER', 'THR', 'TRP', 'TYR', 'VAL']

n_nucl = 4

shutil.copyfile('fiber-167-4_clean.pdb', 'dna.pdb')

histones = helper_functions.parse_pdb('histone-all.pdb')
n_atoms = len(histones.index)

first_atom_per_histone = [0]
last_atom_per_histone = []
for i in range(1, n_atoms):
    c = histones.loc[i, 'chainID']
    prev_c = histones.loc[i - 1, 'chainID']
    if (c == 'A') and (prev_c == 'H'):
        first_atom_per_histone.append(i)
        last_atom_per_histone.append(i - 1)
last_atom_per_histone.append(n_atoms - 1)

assert len(first_atom_per_histone) == n_nucl
assert len(last_atom_per_histone) == n_nucl

for i in range(n_nucl):
    start_id = first_atom_per_histone[i]
    end_id = last_atom_per_histone[i]
    histone_i = histones.loc[start_id:end_id].copy()
    histone_i['serial'] = list(range(len(histone_i.index)))
    helper_functions.write_pdb(histone_i, f'histone_{i}.pdb')

