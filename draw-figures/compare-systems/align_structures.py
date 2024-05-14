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

histone_core = list(range(43, 135, 3)) + list(range(159, 237, 3)) + list(range(257, 352, 3)) + list(range(400, 487, 3)) + list(range(530, 622, 3)) + list(range(646, 724, 3)) + list(range(744, 839, 3)) + list(range(887, 974, 3))
histone_core = np.array(histone_core)
tetra_nucl_167 = mdtraj.load_pdb(pdb_167)
tetra_nucl_172 = mdtraj.load_pdb(pdb_172)
tetra_nucl_177 = mdtraj.load_pdb(pdb_177)

# align
tetra_nucl_172.superpose(tetra_nucl_167, frame=0, atom_indices=histone_core, 
                         ref_atom_indices=histone_core)
tetra_nucl_177.superpose(tetra_nucl_167, frame=0, atom_indices=histone_core, 
                         ref_atom_indices=histone_core)

# save aligned tetra-nucleosomes
tetra_nucl_167.save_pdb(f'{output_dir}/NRL_167_aligned_folded_structure1.pdb')
tetra_nucl_172.save_pdb(f'{output_dir}/NRL_172_aligned_folded_structure1.pdb')
tetra_nucl_177.save_pdb(f'{output_dir}/NRL_177_aligned_folded_structure1.pdb')

# save aligned DNA
# to be continued

