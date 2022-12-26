#########################################################################
# Author: Charlse.Zhang
# Created Time: Thu 03 Apr 2014 05:47:21 PM CDT
# File Name: genConf.sh
# Description: 
#########################################################################
#!/bin/bash

# modified the pdb2cg_dna.py from  ${spn2cdir}/utils to output residue id

export pyPath="/Users/xl23/bin/anaconda2/bin"

PDBname=fiber-${nrl}-${num_copy}_clean.pdb

$pyPath/python ../pdb2cg_dna.py $PDBname

${spn2cdir}/utils/replace_atoms.sh ../../buildDna/bdna_curv_conf_4list.in dna_conf.in dna_from_xtal_template.in

python replace_atom_type.py ../../buildDna/bdna_curv_conf_4list.in dna_from_xtal_template.in dna_from_xtal.in

