# Code for cleaning different data


def write_rot_txt(rot_x, rot_y, rot_z, output_file):
    output_file.write('by rotation x {rot_x}\n'.format(rot_x=rot_x))
    output_file.write('by rotation y {rot_y}\n'.format(rot_y=rot_y))
    output_file.write('by rotation z {rot_z}\n'.format(rot_z=rot_z))


def concat_fake_DNA_with_1zbb_part(input_file_1, input_file_2, output_file, length_fake_DNA):
    # We will remove the residues of fake DNA with index in this list
    remove_resid_idx_list = [length_fake_DNA, length_fake_DNA + 1]
    # we remove the residues with index in remove_resid_idx_list because these residues are used for alignment
    input_file_1_lines = input_file_1.readlines()
    for each_line in input_file_1_lines:
        if each_line[0:4] == 'ATOM' and (int(each_line[22:26]) not in remove_resid_idx_list):
            output_file.write(each_line)
        else:  # leave out the residues of fake DNA with index in remove_resid_idx_list
            pass

    input_file_2_lines = input_file_2.readlines()
    for each_line in input_file_2_lines:
        # Only keep chain I and chain J (chain I and J are DNA chains)
        if each_line[0:4] == 'ATOM' and each_line[21] in ['I', 'J']:
            output_file.write(each_line)
        else:  # leave out histones
            pass


# code for cleaning single nucleosome 
def clean_single_block(input_file, output_file):
    input_file_list = []
    input_file_lines = input_file.readlines()
    # Read the input file and save each line with atom information into a list
    for each_line in input_file_lines:
        if each_line[0:4] == 'ATOM':
            input_file_list.append(each_line)
    N_atom = len(input_file_list)

    # clean chain ID
    for i in range(N_atom):
        if input_file_list[i][21] == 'A':  # Replace chain ID A to I
            input_file_list[i] = input_file_list[i][:21] + 'I' + input_file_list[i][22:]
        if input_file_list[i][21] == 'B':  # Replace chain ID B to J
            input_file_list[i] = input_file_list[i][:21] + 'J' + input_file_list[i][22:]

    atom_index = 1
    residue_index = 1

    # clean chain I
    for i in range(N_atom):
        if input_file_list[i][21] == 'I':
            input_file_list[i] = input_file_list[i][:6] + str(atom_index).rjust(5) + input_file_list[i][11:]  # clean atom serial number
            atom_index += 1

            original_residue_index = int(input_file_list[i][22:26])
            input_file_list[i] = input_file_list[i][:22] + str(residue_index).rjust(4) + input_file_list[i][26:]  # clean residue sequence number

            if i < N_atom - 1 and ((input_file_list[i + 1][21] != input_file_list[i][21]) or (int(input_file_list[i + 1][22:26]) != original_residue_index)):
                residue_index += 1

    # clean chain J
    # we have to divide chain J into two parts
    J_part_1_index_list = []
    J_part_2_index_list = []

    for i in range(N_atom):
        if input_file_list[i][21] == 'J':
            J_part_1_index_list.append(i)
        if input_file_list[i][21] == 'J' and input_file_list[i + 1][21] != 'J':
            break

    for i in range(J_part_1_index_list[-1] + 1, N_atom):
        if input_file_list[i][21] == 'J':
            J_part_2_index_list.append(i)
        if i < N_atom - 1 and input_file_list[i][21] == 'J' and input_file_list[i + 1][21] != 'J':
            break

    # first label the atoms in J_part_2_index_list, then label the atoms in J_part_1_index_list

    for i in J_part_2_index_list:
        input_file_list[i] = input_file_list[i][:6] + str(atom_index).rjust(5) + input_file_list[i][11:]
        atom_index += 1

        original_residue_index = int(input_file_list[i][22:26])
        input_file_list[i] = input_file_list[i][:22] + str(residue_index).rjust(4) + input_file_list[i][26:]

        if i < J_part_2_index_list[-1] and ((input_file_list[i + 1][21] != input_file_list[i][21]) or (int(input_file_list[i + 1][22:26]) != original_residue_index)):
            residue_index += 1

    residue_index += 1  # because the last atom in chain J part 2 is at the end of the file

    for i in J_part_1_index_list:
        input_file_list[i] = input_file_list[i][:6] + str(atom_index).rjust(5) + input_file_list[i][11:]
        atom_index += 1

        original_residue_index = int(input_file_list[i][22:26])
        input_file_list[i] = input_file_list[i][:22] + str(residue_index).rjust(4) + input_file_list[i][26:]

        if i < J_part_1_index_list[-1] and ((input_file_list[i + 1][21] != input_file_list[i][21]) or (int(input_file_list[i + 1][22:26]) != original_residue_index)):
            residue_index += 1

    # write the output according to atom index order
    dict_1 = {}
    for i in range(N_atom):
        # the key of dict_1 is atom serial number
        dict_1[int(input_file_list[i][6:11])] = input_file_list[i]

    # add TER to the end of the chain, and at the ending of the file, there is no TER, but there is an END
    for i in range(1, N_atom + 1):
        output_file.write(dict_1[i])
        if i < N_atom:
            if dict_1[i][21] != dict_1[i + 1][21]:
                output_file.write('TER\n')
    output_file.write('END')


def concat_fiber_with_block(input_file_1, input_file_2, output_file, nrl):
    # input_file_1 is the aligned old chromatin fiber
    # input_file_2 is the aligned block
    # output_file is the new chromatin fiber (one more nucleosome than the old one)
    input_file_2_lines = input_file_2.readlines()
    for each_line in input_file_2_lines:
        if each_line[0:4] == 'ATOM':
            if int(each_line[22:26]) == nrl + 1 or int(each_line[22:26]) == nrl + 2:
                pass
            else:
                output_file.write(each_line)

    input_file_1_lines = input_file_1.readlines()
    for each_line in input_file_1_lines:
        if each_line[0:4] == 'ATOM':
            output_file.write(each_line)
        if each_line[0:3] == 'TER':
            output_file.write(each_line)


def change_chain_id(input_file, output_file):
    input_file_lines = input_file.readlines()
    input_file_line_list = []

    for each_line in input_file_lines:
        if each_line[0:4] == 'ATOM':
            input_file_line_list.append(each_line)

    chain_id_list = [['A', 'B'], ['C', 'D'], ['E', 'F'], ['G', 'H'], ['I', 'J'], ['K', 'L'], ['M', 'N'], ['O', 'P'], ['Q', 'R'], ['S', 'T'], ['U', 'V'], ['W', 'X'],
                     ['Y', 'Z'], ['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h'], ['i', 'j'], ['k', 'l'], ['m', 'n'], ['o', 'p'], ['q', 'r'], ['s', 't'], ['u', 'v'],
                     ['w', 'x'], ['y', 'z'], ['0', '1'], ['2', '3'], ['4', '5'], ['6', '7'], ['8', '9']]  # 62 letters or numbers in all for chain id
    # chain A and I will be given chain id as chain_id_list[:][0]
    # chain B and J will be given chain id as chain_id_list[:][1]

    nucl_index = 0

    for i in range(len(input_file_line_list)):
        if i >= 1:
            if input_file_line_list[i - 1][21] == 'J' and input_file_line_list[i][21] == 'I':
                nucl_index += 1
        if input_file_line_list[i][21] == 'I':
            output_file.write(
                input_file_line_list[i][:21] + chain_id_list[nucl_index][0] + input_file_line_list[i][22:])
        if input_file_line_list[i][21] == 'J':
            output_file.write(
                input_file_line_list[i][:21] + chain_id_list[nucl_index][1] + input_file_line_list[i][22:])


def cut_tail(input_file, output_file, nrl):
    input_file_lines = input_file.readlines()
    for each_line in input_file_lines:
        if ((each_line[21] == 'A' and int(each_line[22:26]) <= nrl - 146) or (each_line[21] == 'B' and int(each_line[22:26]) >= nrl + 149)):
            pass
        else:
            output_file.write(each_line)


def change_atom_index(input_file, output_file, num_copy):
    input_list_1 = []  # save chain A, C, E...
    input_list_2 = [[] for i in range(num_copy)]  # save chain B, D, F...
    chain_I_list = ['A', 'C', 'E', 'G', 'I', 'K', 'M', 'O', 'Q', 'S', 'U', 'W', 'Y', 'a',
                    'c', 'e', 'g', 'i', 'k', 'm', 'o', 'q', 's', 'u', 'w', 'y', '0', '2', '4', '6', '8']
    chain_J_list = ['B', 'D', 'F', 'H', 'J', 'L', 'N', 'P', 'R', 'T', 'V', 'X', 'Z', 'b',
                    'd', 'f', 'h', 'j', 'l', 'n', 'p', 'r', 't', 'v', 'x', 'z', '1', '3', '5', '7', '9']

    # Rearrange the input
    input_file_lines = input_file.readlines()

    for each_line in input_file_lines:
        if each_line[0:4] == 'ATOM':
            chain_id = each_line[21]
            if chain_id in chain_I_list:
                input_list_1.append(each_line)
            if chain_id in chain_J_list:
                input_list_2[chain_J_list.index(chain_id)].append(each_line)

    # Write the output
    atom_index = 1
    residue_index = 1
    old_residue_index = int(input_list_1[0][22:26])
    old_chain_id = input_list_1[0][21]

    # 7-11 for atom serial number (rjust)
    # 22 for chain ID
    # 23-26 for residue sequence number (rjust)

    # First write the output for chain I
    for each in input_list_1:
        new_residue_index = int(each[22:26])
        new_chain_id = each[21]
        if new_residue_index != old_residue_index or old_chain_id != new_chain_id:
            residue_index += 1
        each = each[:6] + str(atom_index).rjust(5) + each[11:21] + 'I' + str(residue_index).rjust(4) + each[26:]
        output_file.write(each)
        atom_index += 1
        old_residue_index = new_residue_index
        old_chain_id = new_chain_id
    output_file.write('TER\n')

    # Then write the output for chain J
    input_list_2.reverse()
    input_list_2_ordered = []
    for each1 in input_list_2:
        for each2 in each1:
            input_list_2_ordered.append(each2)

    for each in input_list_2_ordered:
        new_residue_index = int(each[22:26])
        new_chain_id = each[21]
        if new_residue_index != old_residue_index or old_chain_id != new_chain_id:
            residue_index += 1
        each = each[:6] + str(atom_index).rjust(5) + each[11:21] + 'J' + str(residue_index).rjust(4) + each[26:]
        output_file.write(each)
        atom_index += 1
        old_residue_index = new_residue_index
        old_chain_id = new_chain_id
    output_file.write('END')


def produce_tcl_file(nrl, num_copy, output_file):
    output_file.write('mol new 1kx5.pdb\n')
    output_file.write(
        'mol new fiber-%d-%d_chain_id_changed.pdb\n' % (nrl, num_copy))
    output_file.write('\n')
    output_file.write('set length %d\n' % (nrl))
    output_file.write('set all0 [atomselect 0 "all"]\n')
    output_file.write(
        'set sel0 [atomselect 0 "chain I and name P and resid >= -72 and resid <= 72"]\n')
    output_file.write('set start [expr %d - 144]\n' % (nrl))
    output_file.write('set end [expr %d]\n' % (nrl))
    output_file.write('\n')
    chain_I_list = ['A', 'C', 'E', 'G', 'I', 'K', 'M', 'O', 'Q', 'S', 'U', 'W', 'Y', 'a',
                    'c', 'e', 'g', 'i', 'k', 'm', 'o', 'q', 's', 'u', 'w', 'y', '0', '2', '4', '6', '8']
    chain_J_list = ['B', 'D', 'F', 'H', 'J', 'L', 'N', 'P', 'R', 'T', 'V', 'X', 'Z', 'b',
                    'd', 'f', 'h', 'j', 'l', 'n', 'p', 'r', 't', 'v', 'x', 'z', '1', '3', '5', '7', '9']
    for i in range(1, num_copy + 1):
        output_file.write('set sel%d [atomselect 1 "chain %s and name P and resid >= ${start} and resid <= ${end}"]\n' % (i, chain_I_list[i - 1]))
        output_file.write('set m%d [measure fit $sel0 $sel%d]\n' % (i, i))
        output_file.write('$all0 move $m%d\n' % (i))
        output_file.write(
            'set histone%d [atomselect 0 "not chain I and not chain J"]\n' % (i))
        output_file.write('$histone%d writepdb histone-%d.pdb\n' % (i, i))
        output_file.write('measure rmsd $sel0 $sel%d\n' % (i))
        output_file.write('\n')
    output_file.write('exit')
