import numpy as np
import mdtraj
import sys
import os

"""
Align the pdb structures for better visualization.
"""

output_dir = 'examplary-structures/aligned-structures'
os.makedirs(output_dir, exist_ok=True)

pdb_167 = 'examplary-structures/NRL=167_folded_structure1.pdb'
pdb_172 = 'examplary-structures/NRL=172_folded_structure1.pdb'
pdb_177 = 'examplary-structures/NRL=177_folded_structure1.pdb'

histone_core1_indices = list(range(43, 135, 3)) + list(range(159, 237, 3)) + list(range(257, 352, 3)) + list(range(400, 487, 3)) + list(range(530, 622, 3)) + list(range(646, 724, 3)) + list(range(744, 839, 3)) + list(range(887, 974, 3))
histone_core1_indices = np.array(histone_core1_indices)
tetra_nucl_167 = mdtraj.load_pdb(pdb_167)
tetra_nucl_172 = mdtraj.load_pdb(pdb_172)
tetra_nucl_177 = mdtraj.load_pdb(pdb_177)

# align
tetra_nucl_172.superpose(tetra_nucl_167, frame=0, 
                         atom_indices=histone_core1_indices, 
                         ref_atom_indices=histone_core1_indices)
tetra_nucl_177.superpose(tetra_nucl_167, frame=0, 
                         atom_indices=histone_core1_indices, 
                         ref_atom_indices=histone_core1_indices)

# save aligned tetra-nucleosomes
tetra_nucl_167.save_pdb(f'{output_dir}/NRL_167_aligned_folded_structure1.pdb')
tetra_nucl_172.save_pdb(f'{output_dir}/NRL_172_aligned_folded_structure1.pdb')
tetra_nucl_177.save_pdb(f'{output_dir}/NRL_177_aligned_folded_structure1.pdb')

# save aligned DNA
tetra_nucl_167_dna_indices = np.arange(4 * 974, tetra_nucl_167.n_atoms)
tetra_nucl_172_dna_indices = np.arange(4 * 974, tetra_nucl_172.n_atoms)
tetra_nucl_177_dna_indices = np.arange(4 * 974, tetra_nucl_177.n_atoms)
tetra_nucl_167_dna = tetra_nucl_167.atom_slice(tetra_nucl_167_dna_indices)
tetra_nucl_172_dna = tetra_nucl_172.atom_slice(tetra_nucl_172_dna_indices)
tetra_nucl_177_dna = tetra_nucl_177.atom_slice(tetra_nucl_177_dna_indices)
tetra_nucl_167_dna.save_pdb(f'{output_dir}/NRL_167_aligned_folded_structure1_dna.pdb')
tetra_nucl_172_dna.save_pdb(f'{output_dir}/NRL_172_aligned_folded_structure1_dna.pdb')
tetra_nucl_177_dna.save_pdb(f'{output_dir}/NRL_177_aligned_folded_structure1_dna.pdb')

# save the aligned DNA of the first 3 nucleosomes
n_tetra_nucl_167_bp = 3 * 167 + 147
n_tetra_nucl_172_bp = 3 * 172 + 147
n_tetra_nucl_177_bp = 3 * 177 + 147
n_tri_nucl_167_bp = 2 * 167 + 147
n_tri_nucl_172_bp = 2 * 172 + 147
n_tri_nucl_177_bp = 2 * 177 + 147
selected_indices1 = np.array(list(range(3 * n_tri_nucl_167_bp - 1)) + list(range(6 * n_tetra_nucl_167_bp - 2 - 3 * n_tri_nucl_167_bp, 6 * n_tetra_nucl_167_bp - 2)))
tri_nucl_167_dna = tetra_nucl_167_dna.atom_slice(selected_indices1)
selected_indices2 = np.array(list(range(3 * n_tri_nucl_172_bp - 1)) + list(range(6 * n_tetra_nucl_172_bp - 2 - 3 * n_tri_nucl_172_bp, 6 * n_tetra_nucl_172_bp - 2)))
tri_nucl_172_dna = tetra_nucl_172_dna.atom_slice(selected_indices2)
selected_indices3 = np.array(list(range(3 * n_tri_nucl_177_bp - 1)) + list(range(6 * n_tetra_nucl_177_bp - 2 - 3 * n_tri_nucl_177_bp, 6 * n_tetra_nucl_177_bp - 2)))
tri_nucl_177_dna = tetra_nucl_177_dna.atom_slice(selected_indices3)
tri_nucl_167_dna.save_pdb(f'{output_dir}/NRL_167_aligned_folded_structure1_tri_nucl_dna.pdb')
tri_nucl_172_dna.save_pdb(f'{output_dir}/NRL_172_aligned_folded_structure1_tri_nucl_dna.pdb')
tri_nucl_177_dna.save_pdb(f'{output_dir}/NRL_177_aligned_folded_structure1_tri_nucl_dna.pdb')

