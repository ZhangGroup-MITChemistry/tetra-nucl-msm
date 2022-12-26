#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Dec  1 20:27:23 2021
# File Name: cmd.add_ions.sh
# Description: Add ions based on the current data.prot_dna2 and 
# the desired ionic concentration
#########################################################################
#!/bin/bash

pyPath=$HOME/bin/anaconda2/bin

input_data_file=$1
output_data_file=$2


num_NA=104
num_MG=650
num_CL=0
num_existing_atoms_orig=15564

python getSection.py $input_data_file tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt

num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=68
resid_sta=10384


#python place_ions_0mM_NaCl_5mM_MgCl2.py $num_NA NA_ION.txt -100.0 -300.0 -200.0 59.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0
python place_ions_0mM_NaCl_5mM_MgCl2.py $num_NA NA_ION.txt -100.0 -300.0 -199.0 5.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0

#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_NA" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,35,1.0,$6+5.5,$7+0,$8+0}' > NA_ION.txt

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))


#python place_ions_0mM_NaCl_5mM_MgCl2.py $num_MG MG_ION.txt -96.0 -296.0 -196.0 59.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0
python place_ions_0mM_NaCl_5mM_MgCl2.py $num_MG MG_ION.txt -99.0 -299.0 -198.0 5.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_MG" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,36,2.0,$6+0,$7+5.5,$8+0}' > MG_ION.txt

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

#python place_ions_0mM_NaCl_5mM_MgCl2.py $num_CL CL_ION.txt -93.0 -293.0 -193.0 59.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0
python place_ions_0mM_NaCl_5mM_MgCl2.py $num_CL CL_ION.txt -98.0 -298.0 -197.0 5.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_CL" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,37,-1.0,$6+0,$7+5.5,$8+1}' > CL_ION.txt

num_existing_atoms_updated=$((num_existing_atoms+num_CL))


# Merge them together
cat tmp.txt NA_ION.txt MG_ION.txt CL_ION.txt > modified_ATOM.txt


# Update the original data.prot_dna2 file to include the explicit ions
python updateData_Atoms.py $input_data_file modified_ATOM.txt $output_data_file

# Replace the total number of atoms
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/37 atom types/g; s/-1000.000000 1000.000000 xlo xhi/-100.0 500.0 xlo xhi/g; s/-1000.000000 1000.000000 ylo yhi/-300.0 300.0 ylo yhi/g; s/-1000.000000 1000.000000 zlo zhi/-200.0 400.0 zlo zhi/g; s/34   99.1326/34   99.1326\n  35   22.9898\n  36   24.3050\n  37   35.4530/g" $output_data_file


# Convert to the 4vmd file
$pyPath/python convert_lammps_datafile.py $output_data_file $output_data_file.4vmd
