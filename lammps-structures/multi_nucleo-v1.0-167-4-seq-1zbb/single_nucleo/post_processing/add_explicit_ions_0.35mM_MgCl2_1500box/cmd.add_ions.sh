#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Dec  1 20:27:23 2021
# File Name: cmd.add_ions.sh
# Description: Add ions based on the current data.prot_dna and 
# the desired ionic concentration
#########################################################################
#!/bin/bash

num_NA=0
num_MG=1062
num_CL=1422
num_existing_atoms_orig=7782

python getSection.py data.prot_dna tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt

num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=34
resid_sta=5192


python place_ions.py $num_NA NA_ION.txt -749.0 -749.0 -749.0 70.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_NA" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,35,1.0,$6+5.5,$7+0,$8+0}' > NA_ION.txt

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))


python place_ions.py $num_MG MG_ION.txt -739.0 -739.0 -739.0 70.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_MG" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,36,2.0,$6+0,$7+5.5,$8+0}' > MG_ION.txt

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

python place_ions.py $num_CL CL_ION.txt -729.0 -729.0 -729.0 70.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_CL" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,37,-1.0,$6+0,$7+5.5,$8+1}' > CL_ION.txt


# Merge them together
cat tmp.txt NA_ION.txt MG_ION.txt CL_ION.txt > modified_ATOM.txt


python getSection.py data.prot_dna.4vmd tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt


num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=34
resid_sta=5192

python place_ions_4vmd.py $num_NA NA_ION_4vmd.txt -749.0 -749.0 -749.0 70.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_NA" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,35,1.0,$5+5.5,$6+0,$7+0}' > NA_ION_4vmd.txt

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))

python place_ions_4vmd.py $num_MG MG_ION_4vmd.txt -739.0 -739.0 -739.0 70.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_MG" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,36,2.0,$5+0,$6+5.5,$7+0}' > MG_ION_4vmd.txt

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

python place_ions_4vmd.py $num_CL CL_ION_4vmd.txt -729.0 -729.0 -729.0 70.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_CL" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,37,-1.0,$5+0,$6+0,$7+5.5}' > CL_ION_4vmd.txt

num_existing_atoms_updated=$((num_existing_atoms+num_CL))

# Merge them together
cat tmp.txt NA_ION_4vmd.txt MG_ION_4vmd.txt CL_ION_4vmd.txt > modified_ATOM_4vmd.txt


# Update the original data.prot_dna file to include the explicit ions
python updateData_Atoms.py data.prot_dna modified_ATOM.txt data.prot_dna_ions
python updateData_Atoms.py data.prot_dna.4vmd modified_ATOM_4vmd.txt data.prot_dna_ions.4vmd

# Replace the total number of atoms
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/38 atom types/g; s/-10000.000000 10000.000000 xlo xhi/-750.0 750.0 xlo xhi/g; s/-10000.000000 10000.000000 ylo yhi/-750.0 750.0 ylo yhi/g; s/-10000.000000 10000.000000 zlo zhi/-750.0 750.0 zlo zhi/g; s/34   99.1326/34   99.1326\n  35   22.9898\n  36   24.3050\n  37   35.4530\n  38   49.4000/g" data.prot_dna_ions
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/38 atom types/g; s/-10000.000000 10000.000000 xlo xhi/-750.0 750.0 xlo xhi/g; s/-10000.000000 10000.000000 ylo yhi/-750.0 750.0 ylo yhi/g; s/-10000.000000 10000.000000 zlo zhi/-750.0 750.0 zlo zhi/g; s/34   99.1326/34   99.1326\n  35   22.9898\n  36   24.3050\n  37   35.4530\n  38   49.4000/g" data.prot_dna_ions.4vmd

