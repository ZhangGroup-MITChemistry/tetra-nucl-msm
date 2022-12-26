# Use this python script to build the chromatin fiber with different NRL and number of nucleosomes

import os
import tools
import sys

# Notice there are some limitaions for this program:
# (1) We can only build chromatin fiber with fixed NRL
# (2) Cannot build chromatin fiber with super long linker DNA or super short linker DNA
# (3) If the number of nucleosomes is larger than 31 (i.e. num_copy > 31), this program may not be able to plug in histones, since no enough letter or number for chain ID

# 0-Define some parameters (all the input parameters are in this section)
nrl = int(sys.argv[1])  
seq_fake_DNA = sys.argv[2]  
num_copy = int(sys.argv[3])  
rot_x, rot_y, rot_z = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]) 

# 1-build the "fake" DNA chain with given sequence
length_fake_DNA = nrl - 165
if len(seq_fake_DNA) != length_fake_DNA:  # check if the length of fake DNA is correct
    print('The length of fake DNA is not correct! (i.e. not consistent with the NRL)')

os.system('fiber -b -seq=%s %d-bp.pdb' % (seq_fake_DNA, length_fake_DNA))

# 2-build the first monomer
# 2.1-Align 1zbb_part-167.pdb and fake DNA
os.system('find_pair 1zbb_part-167.pdb stdout')
os.system('frame_mol -1 ref_frames.dat 1zbb_part-167.pdb 1zbb_part-167_aligned.pdb')
os.system('find_pair {seq_length}-bp.pdb {seq_length}-bp.bps'.format(seq_length=length_fake_DNA))
os.system('frame_mol -{seq_length} ref_frames.dat {seq_length}-bp.pdb {seq_length}-bp_aligned.pdb'.format(seq_length=length_fake_DNA))
output_file = open('rot_xyz.txt', 'w')
tools.write_rot_txt(rot_x, rot_y, rot_z, output_file)
output_file.close()
os.system('rotate_mol -r=rot_xyz.txt {seq_length}-bp_aligned.pdb {seq_length}-bp_aligned_rotated.pdb'.format(seq_length=length_fake_DNA))

# Then we concatenate the fake DNA with 1zbb_part-167.pdb
input_file_1 = open('%d-bp_aligned_rotated.pdb' % (length_fake_DNA), 'r')
input_file_2 = open('1zbb_part-167_aligned.pdb', 'r')
# Each monomer has (nrl + 1) base pairs
output_file = open('block-%d.pdb' % (nrl + 1), 'w')

tools.concat_fake_DNA_with_1zbb_part(input_file_1, input_file_2, output_file, length_fake_DNA)

input_file_1.close()
input_file_2.close()
output_file.close()

# 2.2-Clean the monomer structure file (this part is super dirty!!!)
# So our problem is, the fake DNA has chain ID A or B, while DNA in 1zbb_part-167.pdb has chain ID I or J
# We have to rearrange the chain ID and residue index of each atom
# block-{block_length}-clean.pdb is our building block
input_file = open('block-%d.pdb' % (nrl + 1), 'r')
output_file = open('block-%d-clean.pdb' % (nrl + 1), 'w')
# use this function to do the job
tools.clean_single_block(input_file, output_file)
input_file.close()
output_file.close()

# 3-build the chromatin fiber with num_copy nucleosomes
# fiber-{length}-n.pdb is the chromatin fiber with NRL = length and n nucleosomes
# for block-{block_length}-clean.pdb, block_length = NRL + 1, since one extra bp for alignment and overlap
os.system('cp block-%d-clean.pdb fiber-%d-1.pdb' % (nrl + 1, nrl))

for n in range(1, num_copy):
    os.system('find_pair fiber-%d-%d.pdb stdout' % (nrl, n))
    os.system('frame_mol -1 ref_frames.dat fiber-%d-%d.pdb fiber-%d-%d_aligned.pdb' % (nrl, n, nrl, n))
    os.system('find_pair block-%d-clean.pdb block-%d-clean.bps' % (nrl + 1, nrl + 1))
    os.system('frame_mol -%d ref_frames.dat block-%d-clean.pdb block-%d-clean_aligned.pdb' % (nrl, nrl + 1, nrl + 1))
    # note here we use frame_mol -nrl, instead of -(nrl + 1), since there is one bp that 3dna cannot recognize
    # that is to say, 3dna recognizes that block-(nrl+1)-clean.pdb has only nrl bps, not (nrl+1) bps
    # check block-(nrl+1)-clean.bps and you can find that one bp is not recognized
    # actually, for 1zbb_part-167.pdb, 3dna can only recognize 166 bp of it
    # you can run "find_pair 1zbb_part-167.pdb 1zbb_part-167.bps" and check 1zbb_part-167.bps, and you will see only 166 bps are recognized

    # Then we concatenate the block with the n-mer fiber to get (n+1)-mer fiber
    input_file_1 = open('fiber-%d-%d_aligned.pdb' % (nrl, n), 'r')
    input_file_2 = open('block-%d-clean_aligned.pdb' % (nrl + 1), 'r')
    output_file = open('fiber-%d-%d.pdb' % (nrl, n + 1), 'w')
    tools.concat_fiber_with_block(input_file_1, input_file_2, output_file, nrl)
    input_file_1.close()
    input_file_2.close()
    output_file.close()
    os.system('rm -f fiber-%d-%d.pdb' % (nrl, n))
    os.system('rm -f fiber-%d-%d_aligned.pdb' % (nrl, n))

# 4-change chain id, cut the tail and re-index the chromatin fiber

# change chain id, so that we can align histones
input_file = open('fiber-%d-%d.pdb' % (nrl, num_copy), 'r')
output_file = open('fiber-%d-%d_chain_id_changed_with_tail.pdb' % (nrl, num_copy), 'w')
tools.change_chain_id(input_file, output_file)
input_file.close()
output_file.close()
os.system('rm -f fiber-%d-%d.pdb' % (nrl, num_copy))

# cut the tail in chain A and chain B
input_file = open('fiber-%d-%d_chain_id_changed_with_tail.pdb' % (nrl, num_copy), 'r')
output_file = open('fiber-%d-%d_chain_id_changed.pdb' % (nrl, num_copy), 'w')
tools.cut_tail(input_file, output_file, nrl)
input_file.close()
output_file.close()

# re-index the fiber (i.e. change the index of atoms)
input_file = open('fiber-%d-%d_chain_id_changed.pdb' % (nrl, num_copy), 'r')
output_file = open('fiber-%d-%d_clean.pdb' % (nrl, num_copy), 'w')
tools.change_atom_index(input_file, output_file, num_copy)
input_file.close()
output_file.close()

# 5-plug in histones
# Produce the .tcl file automatically to further plug in histones with VMD
output_file = open('align_histone.tcl', 'w')
tools.produce_tcl_file(nrl, num_copy, output_file)
output_file.close()
