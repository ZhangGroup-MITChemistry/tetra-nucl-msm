#########################################################################
# Author: Xingcheng Lin
# Created Time: Fri Mar 18 00:00:02 2022
# File Name: cmd.add_ions_to_converted_data_file.sh 
# Description: This script will add ions to the converted lammps data input file
# for the next stage of simulations
#########################################################################
#!/bin/bash

export startingRef1=2.5
export startingRef2=50.0

startID=$1
endID=$2

for ((i=0; i<=2; i++));
do
    for ((j=0; j<=4; j+=1));
    do
        ref1=$(bc<<<"$startingRef1+2.5*($i)" | awk '{printf "%1.1f\n", $0}')
        ref2=$(bc<<<"$startingRef2+50.0*($j)" | awk '{printf "%1.1f\n", $0}')
        
        # Add ions to the converted lammps data file
#        mkdir -p starting_structure_for2copies/15mM_NaCl_5mM_MgCl2/
#        bash cmd.add_ions_15mM_NaCl_5mM_MgCl2.sh starting_structure_for2copies/data.prot_dna_${ref1}_${ref2} starting_structure_for2copies/15mM_NaCl_5mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
#        # Convert to a dcd file
#        vmd -e convert_to_dcd.tcl -args starting_structure_for2copies/15mM_NaCl_5mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
#       
#        mkdir -p starting_structure_for2copies/1.5mM_NaCl_5mM_MgCl2/
#        bash cmd.add_ions_1.5mM_NaCl_5mM_MgCl2.sh starting_structure_for2copies/data.prot_dna_${ref1}_${ref2} starting_structure_for2copies/1.5mM_NaCl_5mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
#        # Convert to a dcd file
#        vmd -e convert_to_dcd.tcl -args starting_structure_for2copies/1.5mM_NaCl_5mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
#
#
#        mkdir -p starting_structure_for2copies/150mM_NaCl_0mM_MgCl2/
#        bash cmd.add_ions_150mM_NaCl_0mM_MgCl2.sh starting_structure_for2copies/data.prot_dna_${ref1}_${ref2} starting_structure_for2copies/150mM_NaCl_0mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
#        # Convert to a dcd file
#        vmd -e convert_to_dcd.tcl -args starting_structure_for2copies/150mM_NaCl_0mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
#
        mkdir -p starting_structure_for2copies/0mM_NaCl_5mM_MgCl2/
        bash cmd.add_ions_0mM_NaCl_5mM_MgCl2.sh starting_structure_for2copies/data.prot_dna_${ref1}_${ref2} starting_structure_for2copies/0mM_NaCl_5mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}
        # Convert to a dcd file
        vmd -e convert_to_dcd.tcl -args starting_structure_for2copies/0mM_NaCl_5mM_MgCl2/data.prot_dna_ions_${ref1}_${ref2}


    done
done

