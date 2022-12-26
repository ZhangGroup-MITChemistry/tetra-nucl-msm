#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Dec  1 20:27:23 2021
# File Name: cmd.add_ions.sh
# Description: Add ions based on the current data.prot_dna2 and 
# the desired ionic concentration
#########################################################################
#!/bin/bash

num_NA=19505
num_MG=650
num_CL=19401
num_existing_atoms_orig=15564

python getSection.py data.prot_dna2 tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt

num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=68
resid_sta=10384


python place_ions.py $num_NA NA_ION.txt -100.0 -300.0 -200.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_NA" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,35,1.0,$6+5.5,$7+0,$8+0}' > NA_ION.txt

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))


python place_ions.py $num_MG MG_ION.txt -96.0 -296.0 -196.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_MG" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,36,2.0,$6+0,$7+5.5,$8+0}' > MG_ION.txt

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

python place_ions.py $num_CL CL_ION.txt -93.0 -293.0 -193.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_CL" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,residsta+NR,37,-1.0,$6+0,$7+5.5,$8+1}' > CL_ION.txt


# Merge them together
cat tmp.txt NA_ION.txt MG_ION.txt CL_ION.txt > modified_ATOM.txt


python getSection.py data.prot_dna2.4vmd tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt


num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=68
resid_sta=10384

python place_ions_4vmd.py $num_NA NA_ION_4vmd.txt -100.0 -300.0 -200.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_NA" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,35,1.0,$5+5.5,$6+0,$7+0}' > NA_ION_4vmd.txt

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))

python place_ions_4vmd.py $num_MG MG_ION_4vmd.txt -96.0 -296.0 -196.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_MG" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,36,2.0,$5+0,$6+5.5,$7+0}' > MG_ION_4vmd.txt

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

python place_ions_4vmd.py $num_CL CL_ION_4vmd.txt -93.0 -293.0 -193.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0
#cat tmp.txt | awk -v chainIDsta="$chainID_sta" -v residsta="$resid_sta" -v numION="$num_CL" -v numExistingAtoms="$num_existing_atoms" '{if(NR<=numION) print " ",numExistingAtoms+NR,chainIDsta+NR,37,-1.0,$5+0,$6+0,$7+5.5}' > CL_ION_4vmd.txt

num_existing_atoms_updated=$((num_existing_atoms+num_CL))

# Merge them together
cat tmp.txt NA_ION_4vmd.txt MG_ION_4vmd.txt CL_ION_4vmd.txt > modified_ATOM_4vmd.txt


# Update the original data.prot_dna2 file to include the explicit ions
python updateData_Atoms.py data.prot_dna2 modified_ATOM.txt data.prot_dna2_ions
python updateData_Atoms.py data.prot_dna2.4vmd modified_ATOM_4vmd.txt data.prot_dna2_ions.4vmd

# Replace the total number of atoms
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/37 atom types/g; s/-1000.000000 1000.000000 xlo xhi/-100.0 500.0 xlo xhi/g; s/-1000.000000 1000.000000 ylo yhi/-300.0 300.0 ylo yhi/g; s/-1000.000000 1000.000000 zlo zhi/-200.0 400.0 zlo zhi/g" data.prot_dna2_ions
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/37 atom types/g; s/-1000.000000 1000.000000 xlo xhi/-100.0 500.0 xlo xhi/g; s/-1000.000000 1000.000000 ylo yhi/-300.0 300.0 ylo yhi/g; s/-1000.000000 1000.000000 zlo zhi/-200.0 400.0 zlo zhi/g" data.prot_dna2_ions.4vmd

