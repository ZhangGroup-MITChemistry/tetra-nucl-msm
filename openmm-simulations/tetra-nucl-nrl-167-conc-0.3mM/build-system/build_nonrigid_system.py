import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj

sys.path.append('/home/gridsan/sliu/Projects/smog-3spn2-openmm')

from OpenSMOG3SPN2.forcefields.parsers import SMOGParser, DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq, write_pdb
from OpenSMOG3SPN2.utils.insert import insert_molecules

"""
Build system xml file without any rigid bodies.
In this system, bonds and angles are all kept, and exclusions are all kept. 
"""

if not os.path.exists('intermediate-files'):
    os.makedirs('intermediate-files')

n_single_nucl = 26
box_a, box_b, box_c = 55.0, 55.0, 55.0
chromatin = SMOG3SPN2Model()

# load tetranucleosome histones
for i in range(4):
    histone_i_parser = SMOGParser.from_atomistic_pdb(f'../tetra-nucl-nrl-167-pdb-files/histone_{i}.pdb', 
                                                     f'intermediate-files/cg_tetra_nucl_histone_{i}.pdb', 
                                                     default_parse=False)
    histone_i_parser.parse_mol(get_native_pairs=False)
    chromatin.append_mol(histone_i_parser)

# load tetranucleosome DNA
with open('tetra_nucl_nrl_167_dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
full_seq1 = seq1 + get_WC_paired_seq(seq1)

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../tetra-nucl-nrl-167-pdb-files/dna.pdb', 
                                               'intermediate-files/cg_tetra_nucl_dna.pdb', new_sequence=full_seq1, 
                                               temp_name='dna1')
chromatin.append_mol(dna_parser)

chromatin.atoms_to_pdb('intermediate-files/cg_tetra_nucl_original_coord.pdb')
n_tetra_nucl_atoms = len(chromatin.atoms.index)
print(f'{n_tetra_nucl_atoms} CG atoms in tetranucleosome. ')

# put tetranucleosome geometric center at the center of the box
tetra_nucl_atoms = chromatin.atoms.copy()
coord = tetra_nucl_atoms[['x', 'y', 'z']].to_numpy()
coord -= np.mean(coord, axis=0)
coord += 10*np.array([box_a, box_b, box_c])/2
tetra_nucl_atoms[['x', 'y', 'z']] = coord
tetra_nucl_atoms.loc[:, 'charge'] = ''
write_pdb(tetra_nucl_atoms, 'cg_tetra_nucl.pdb')

# load single nucleosomes, first load histone, then load DNA
histone_parser = SMOGParser.from_atomistic_pdb('../single-nucl-pdb-files/histone.pdb', 
                                               'intermediate-files/cg_single_nucl_histone.pdb', default_parse=False)
histone_parser.parse_mol(get_native_pairs=False)
with open('single_nucl_dna_seq.txt', 'r') as f:
    seq2 = f.readlines()[0].strip()
full_seq2 = seq2 + get_WC_paired_seq(seq2)
dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../single-nucl-pdb-files/dna.pdb', 
                                               'intermediate-files/cg_single_nucl_dna.pdb', 
                                               new_sequence=full_seq2, temp_name='dna2')
nucl = SMOG3SPN2Model()
nucl.append_mol(histone_parser)
nucl.append_mol(dna_parser)
nucl.atoms_to_pdb('cg_single_nucl.pdb')
n_atoms_per_single_nucl = len(nucl.atoms.index)
print(f'{n_atoms_per_single_nucl} CG atoms in each single nucleosome. ')

# prepare pdb file for the whole system
if not os.path.exists('start.pdb'):
    insert_molecules('cg_single_nucl.pdb', 'start.pdb', n_mol=n_single_nucl, existing_pdb='cg_tetra_nucl.pdb', 
                     box=[box_a, box_b, box_c])
top = app.PDBFile('start.pdb').getTopology()

for i in range(n_single_nucl):
    chromatin.append_mol(histone_parser)
    chromatin.append_mol(dna_parser)

chromatin.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c)
chromatin.add_protein_bonds(force_group=1)
chromatin.add_protein_angles(force_group=2)
chromatin.add_dna_bonds(force_group=5)
chromatin.add_dna_angles(force_group=6)
chromatin.add_dna_stackings(force_group=7)
chromatin.add_dna_dihedrals(force_group=8)
chromatin.add_dna_base_pairs(force_group=9)
chromatin.add_dna_cross_stackings(force_group=10)
chromatin.parse_all_exclusions()
chromatin.add_all_vdwl(force_group=11)
chromatin.add_all_elec(force_group=12)
chromatin.save_system('nonrigid_system.xml')


