#########################################################################
# Author: Xingcheng Lin
# Created Time: Tue Jul  6 15:00:42 2021
# File Name: cmd_create_212mers.sh
# Description: 
#########################################################################
#!/bin/bash

#vmd -e create_2-12mers.tcl -args "100.0"
vmd -e create_2copies.tcl

cat mono1.pdb mono2.pdb > 2copies.pdb

# Delete the extra lines for the next formatting steps
gsed -i '/CRYST/d' ./2copies.pdb
gsed -i '/END/d' ./2copies.pdb

# Use MayChemTools to renumber residues and atoms
ModifyPDBFiles.pl --mode RenumberAtoms 2copies.pdb 
ModifyPDBFiles.pl --mode RenumberResidues 2copiesRenumberAtoms.pdb 


# Convert to dcd file
mdconvert -f -o starting_structure_for2copies/start_2copies.dcd 2copiesRenumberAtomsRenumberResidues.pdb 
