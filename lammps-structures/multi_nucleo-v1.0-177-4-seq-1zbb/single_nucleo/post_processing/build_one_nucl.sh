#########################################################################
# Author: Xingcheng Lin
# Created Time: Fri Mar  8 21:05:31 2019
# File Name: build_one_nucl.sh
# Description: 
#########################################################################
#!/bin/bash

export pyPath="/Users/xl23/bin/anaconda2/bin"
cd dna_list_files/
$pyPath/python dna_list.py
cd ../

cp ../sbm_scripts/single_nucl/ca_angle_list.txt prot_list_files/
cp ../sbm_scripts/single_nucl/ca_dihed_list.txt prot_list_files/
cp ../sbm_scripts/single_nucl/ca_bond_list.txt prot_list_files/
cp ../sbm_scripts/single_nucl/ca_pair_list.txt prot_list_files/
cp ../sbm_scripts/single_nucl/exclusion_list.txt prot_list_files/

$pyPath/python build_one_nucl.py

# Build rigid group
cd groupsDefinition/
$pyPath/python setup_rigidGroups_nucl.py
cd ../


