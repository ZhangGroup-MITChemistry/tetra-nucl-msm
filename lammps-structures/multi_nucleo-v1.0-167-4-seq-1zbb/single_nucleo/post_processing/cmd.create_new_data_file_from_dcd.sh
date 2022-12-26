#########################################################################
# Author: Xingcheng Lin
# Created Time: Fri Mar 18 00:00:02 2022
# File Name: cmd.create_new_data_file_from_dcd.sh
# Description: This script will convert the dcd file extract our simulation
# and map that into the lammps data file as new input
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
        echo $ref1 $ref2

        gsed "s/CV1/${ref1}/g; s/CV2/${ref2}/g" template_convert_to_lammpsdata.tcl > convert_to_lammpsdata.tcl

        # Create a lammps data file from the simulation-extracted dcd file
        vmd -e convert_to_lammpsdata.tcl

        # update the data file for this converted snapshot
        python updateData_Atoms_fromlammpsdata.py data.prot_dna2 data.tmp starting_structure_for2copies/data.prot_dna_${ref1}_${ref2}
    done
done

