import numpy as np

# ----- input parameters -----
tetra_nucl_nrl = 167
n_single_nucl = 26
# ----------------------------

# note for plumed, atom id starts from 1
histone_core_start_id = [44, 160, 258, 401, 531, 647, 745, 888]
histone_core_end_id = [135, 237, 352, 487, 622, 724, 839, 974]
n_atoms_per_histone = 974

n_tetra_nucl_bp = 3*tetra_nucl_nrl + 147
n_tetra_nucl_histone_atoms = n_atoms_per_histone*4
n_tetra_nucl_dna_atoms = 6*n_tetra_nucl_bp - 2
n_tetra_nucl_atoms = n_tetra_nucl_histone_atoms + n_tetra_nucl_dna_atoms

n_dna_atoms_per_single_nucl = 6*147 - 2
n_atoms_per_single_nucl = n_atoms_per_histone + n_dna_atoms_per_single_nucl

# tetra-nucleosome bp id starts from 1
tetra_nucl_ssDNA1_bp_id_to_atom_id = {}
for i in range(n_tetra_nucl_bp):
    if i == 0:
        atoms1 = (np.arange(2) + n_tetra_nucl_histone_atoms + 1).tolist()
    else:
        atoms1 = (np.arange(3) + n_tetra_nucl_histone_atoms + 3*i).tolist()
    tetra_nucl_ssDNA1_bp_id_to_atom_id[i + 1] = atoms1

tetra_nucl_core_ssDNA1_start_id = []
tetra_nucl_core_ssDNA1_end_id = []
tetra_nucl_linker_ssDNA1_start_id = []
tetra_nucl_linker_ssDNA1_end_id = []
for i in range(4):
    b1 = i*tetra_nucl_nrl + 1
    b2 = b1 + 147 - 1
    a1 = min(tetra_nucl_ssDNA1_bp_id_to_atom_id[b1])
    a2 = max(tetra_nucl_ssDNA1_bp_id_to_atom_id[b2])
    tetra_nucl_core_ssDNA1_start_id.append(a1)
    tetra_nucl_core_ssDNA1_end_id.append(a2)
    if i <= 2:
        b3 = b2 + 1
        b4 = b3 + (tetra_nucl_nrl - 147) - 1
        a3 = min(tetra_nucl_ssDNA1_bp_id_to_atom_id[b3])
        a4 = max(tetra_nucl_ssDNA1_bp_id_to_atom_id[b4])
        tetra_nucl_linker_ssDNA1_start_id.append(a3)
        tetra_nucl_linker_ssDNA1_end_id.append(a4)

with open('plumed.txt', 'w') as f:
    f.write('UNITS LENGTH=nm TIME=ps ENERGY=kj/mol\n\n')
    
    # rebuild tetra-nucleosome
    entity0 = []
    for i in range(4):
        a1 = tetra_nucl_core_ssDNA1_start_id[i]
        a2 = tetra_nucl_core_ssDNA1_end_id[i]
        entity0.append(f'{a1}-{a2}')
        a3 = i*n_atoms_per_histone + 1
        a4 = a3 + n_atoms_per_histone - 1
        entity0.append(f'{a3}-{a4}')
        if i <= 2:
            a5 = tetra_nucl_linker_ssDNA1_start_id[i]
            a6 = tetra_nucl_linker_ssDNA1_end_id[i]
            entity0.append(f'{a5}-{a6}')
    a7 = n_tetra_nucl_histone_atoms + 3*n_tetra_nucl_bp
    a8 = n_tetra_nucl_atoms
    entity0.append(f'{a7}-{a8}')
    entity0 = ','.join(entity0)
    f.write(f'WHOLEMOLECULES ENTITY0={entity0}\n\n')
    
    for i in range(4):
        sel = []
        for j in range(len(histone_core_start_id)):
            a1 = histone_core_start_id[j] + i*n_atoms_per_histone
            a2 = histone_core_end_id[j] + i*n_atoms_per_histone
            sel.append(f'{a1}-{a2}:3')
        sel = ','.join(sel)
        f.write(f'nucl{i + 1}: CENTER ATOMS={sel} NOPBC\n')
    for i in range(n_single_nucl):
        sel = []
        for j in range(len(histone_core_start_id)):
            a1 = histone_core_start_id[j] + n_tetra_nucl_atoms + i*n_atoms_per_single_nucl
            a2 = histone_core_end_id[j] + n_tetra_nucl_atoms + i*n_atoms_per_single_nucl
            sel.append(f'{a1}-{a2}:3')
        sel = ','.join(sel)
        f.write(f'nucl{i + 5}: CENTER ATOMS={sel}\n')
    f.write('\n')
    for i in range(3):
        for j in range(i + 1, 4):
            f.write(f'd{i + 1}{j + 1}: DISTANCE ATOMS=nucl{i + 1},nucl{j + 1} NOPBC\n')
    f.write('\n')
    sel = []
    for i in range(3):
        for j in range(i + 1, 4):
            sel.append(f'd{i + 1}{j + 1}')
    sel = ','.join(sel)
    f.write(f'PRINT STRIDE=500 ARG={sel} FILE=COLVARS\n')
    sel = [f'nucl{i + 1}' for i in range(4 + n_single_nucl)]
    sel = ','.join(sel)
    f.write(f'DUMPATOMS STRIDE=500 FILE=nucl.xyz ATOMS={sel}\n')
    f.write(f'FLUSH STRIDE=1000\n\n')

