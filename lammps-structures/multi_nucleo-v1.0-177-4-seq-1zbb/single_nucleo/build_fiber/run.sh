#!/bin/bash

# Run this script to build the chromatin fiber
# You only need to change the parameters in parameter section in this script, and then bash this script

# The files you need: 
# (a) this run.sh script
# (b) main-1.py & main-2.py
# (c) tools.py
# (d) 1kx5.pdb
# (e) 1zbb_part-167.pdb

# running section
# if you do not know what you are doing, please do not change the code below
# Warning: if you receive error "vmd: command not found", then run order "vmd -dispdev text -e align_histone.tcl" in terminal manually

# ——————————————run the orders—————————————————
python main-1.py $nrl $sequence $num_copy $x_rot $y_rot $z_rot
vmd -dispdev text -e align_histone.tcl
python main-2.py $nrl $num_copy
# ——————————————end—————————————————

# The final outputs are fiber-{nrl}-{num_copy}_clean.pdb & histone-all.pdb



