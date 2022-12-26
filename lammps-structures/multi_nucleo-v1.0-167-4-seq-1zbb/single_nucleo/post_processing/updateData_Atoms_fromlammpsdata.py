####################################################################################
# This script will help read in the data.prot_dna file and output a new dara.prot_dna file 
# based on the updated Atom section from the converted lammps data file from vmd topo
# writelammpsdata tool (https://sites.google.com/site/akohlmey/software/topotools/documentation?authuser=0);

# Written by Xingcheng Lin, 03/17/2022
####################################################################################

import math;
import subprocess;
import os;
import math;
import numpy as np;
import sys;

################################################
def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step

def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
###########################################

def updateTop(inputTopFile, editedLammpsFile, outputTopFile):


    infile1 = open(inputTopFile, 'r');
    infile2 = open(editedLammpsFile, 'r');

    outfile = open(outputTopFile, 'w');

    lines1 = [line.rstrip() for line in infile1];
    lines2 = [line.rstrip() for line in infile2];

    length1 = len(lines1);
    length2 = len(lines2);

    
    for i in my_lt_range(0, length1, 1):
       
        line1 = lines1[i].split()

        try:
            print(line1[0])
        except IndexError:
            continue
        else:
            if (line1[0] == "Atoms"):
                idx_atoms = i
            elif (line1 [0] == "Bonds"):
                idx_bonds = i

    for j in my_lt_range(0, length2, 1):
       
        line2 = lines2[j].split()

        try:
            print(line2[0])
        except IndexError:
            continue
        else:
            if (line2[0] == "Atoms"):
                jdx_atoms = j
            elif (line2[0] == "Bonds"):
                jdx_bonds = j

    # Total number of lines for replacing x, y, z coordinates
    num_replacement_lines = idx_bonds - idx_atoms
    # Start copying
    # Copying flag
    flag1 = 0

    for i in my_lt_range(0, length1, 1):
       
        line1 = lines1[i].split()

        if (i <= idx_atoms):
            outfile.write(lines1[i] + "\n")
        elif (i > idx_atoms and i < idx_bonds):

            try:
                print(line1[0])
            except IndexError:
                outfile.write("\n")
                continue
            else:
                delta_idx = i - idx_atoms
                j = jdx_atoms + delta_idx
                line2 = lines2[j].split()

                # If atom id and chain id is the same, it is the same atom, the topo tool messed up the residue ID, so ignoring them
                if (line2[0] == line1[0] and line2[1] == line1[1]):
                    # The 5,6,7 columns of the converted lammps data file is the x, y, z coordinates, while the 6,7,8 columns of the original data file is the x, y, z coordinates
                    outfile.write(line1[0] + " " + line1[1] + " " + line1[2] + " " + line1[3] + " " + line1[4] + " " + line2[4] + " " + line2[5] + " " + line2[6] + "\n")
        else:
            outfile.write(lines1[i] + "\n")
                

############################################################################

if __name__ == "__main__":

    inputTopFile = sys.argv[1]
    editedLammpsFile = sys.argv[2]
    outputTopFile = sys.argv[3]

    updateTop(inputTopFile, editedLammpsFile, outputTopFile);


    print("Love is an endless mystery,")
    print("for it has nothing else to explain it.")

