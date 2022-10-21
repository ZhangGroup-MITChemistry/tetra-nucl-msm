# Notes and instructions about how to run this workflow
# Though the name of some files and folders in multi_nucleo is named for single nucleosome system, actually you can use this to build the multi-nucleosome fiber
# It is only because I modify the code and scripts originally for building single nucleosome, and I do not change the name of the folders and scripts to avoid some path problem
# You need to manually change the path of lmp_tools and some other packages 

# Instructions:
# To build chromatin fiber with given nrl and number of nucleosomes, you only need to do these 4 things:
# (1) First manually change the parameters and paths in this script
# (2) Manually put dnaSeq.txt at single_nucleo/buildDna. Notice that dnaSeq.txt has nrl * (num_copy - 1) + 147 bps
# (3) Manually modify single_nucleo/sbm_scripts/convert_gromacs_to_lammps.py, single_nucleo/post_processing/dna_list_files/dna_list.py, single_nucleo/post_processing/merge_mols.py and single_nucleo/post_processing/groupsDefinition/setup_rigidGroups_nucl.py to add the path for lmp_tools/lib
# (4) Finally run this script (under the path of where this script is)

# A useful tip is, if the structure you build has many overlaps, you may receive some errors from smog. One way to ignore such errors is to add "-warn -1" in file single_nucleo/sbm_scripts/single_nucl/smog2/cmd.smog2.sh (check smog manual for more details)

# ——————————————parameters & paths—————————————————
# nrl is the nucleosome repeat length
# sequence is the sequence of the extended fake DNA, and its length should be NRL - 165
# num_copy is the number of nucleosomes
# x_rot, y_rot and z_rot are the rotation angles of fake DNA relative to 1zbb_part-167.pdb
# if you find overlap for your chromatin, you should consider adjusting x_rot, y_rot and z_rot
nrl=177
sequence=ggcggccagtac
num_copy=4
x_rot=0
y_rot=15
z_rot=50

export smog2dir=/Users/xl23/bin/smog-2.2
spn2cdir="/Users/xl23/GitHub/USER-3SPN2"
export X3DNA="/Users/xl23/bin/x3dna-v2.4"
export PATH="/Users/xl23/bin/x3dna-v2.4/bin:$PATH"
export pyPath="/Users/xl23/bin/anaconda2/bin"
# also assign the path for python2
# ——————————————end—————————————————

# build the fiber structures and put the output structures into corresponding folders
cd single_nucleo/build_fiber/
source run.sh
cp histone-all.pdb ../sbm_scripts/single_nucl/smog2/
cp fiber-${nrl}-${num_copy}_clean.pdb ../post_processing/dnaCoord_from_xtal # Notice that nrl and num_copy is defined in single_nucleo/build_fiber/run.sh

# remove IDR and add ff to histones 
cd ../sbm_scripts/single_nucl/smog2/
source cmd.smog2.sh
cp smog.gro ..
cp smog.top ..
cd ../
source cmd.remove_idr.sh
$pyPath/python gen_lmps_file.py $num_copy

# produce a DNA chain with given sequence
cd ../../buildDna/
source genConf.sh

# post-process (combine DNA and histone ff)
cd ../post_processing/dnaCoord_from_xtal/
source genConf.sh
cd ..
source build_one_nucl.sh

