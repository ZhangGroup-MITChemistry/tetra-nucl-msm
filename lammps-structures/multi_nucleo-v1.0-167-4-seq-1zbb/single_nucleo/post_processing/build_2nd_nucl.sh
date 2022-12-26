#########################################################################
# Author: Xingcheng Lin
# Created Time: Fri Mar  8 21:05:31 2019
# File Name: build_2nd_nucl.sh
# Description: 
#########################################################################
#!/bin/bash

pyPath="/Users/xl23/bin/anaconda2/bin/"

# Replicate DNA for the DNA topology file
cd dna_list_files_replica/
$pyPath/python dna_list_replicate.py
cd ../

# Replicate protein for the protein topology file
cd prot_list_files_replica/
$pyPath/python prot_list_replicate.py
cd ../

$pyPath/python build_2nd_nucl.py

# Build rigid group
cd groupsDefinition/
$pyPath/python setup_rigidGroups_2nucl.py
cd ../


