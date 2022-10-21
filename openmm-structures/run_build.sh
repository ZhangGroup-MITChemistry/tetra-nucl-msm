# build the tetranucleosome structure for openmm simulations

# set some paths and parameters
ca_sbm_3spn2c_openmm=/home/gridsan/sliu/Projects/chromatin-aggregation/CA_SBM_3SPN2C_OPENMM
n_nucl=4
nrl=177
lammps_data_main_dir=/home/gridsan/sliu/Projects/tetra-nucl-msm/lammps-structures/multi_nucleo-v1.0-${nrl}-${n_nucl}-seq-1zbb

# set working directory
mkdir -p nrl-${nrl}-${n_nucl}mer
cd nrl-${nrl}-${n_nucl}mer

# prepare files
# prepare dna and protein pdb files
mkdir -p pdb-files
cd pdb-files
cp ${lammps_data_main_dir}/single_nucleo/build_fiber/fiber-${nrl}-${n_nucl}_clean.pdb .
cp ${lammps_data_main_dir}/single_nucleo/build_fiber/histone-all.pdb .
python ${ca_sbm_3spn2c_openmm}/scripts/separate_histone_dna.py --dna fiber-${nrl}-${n_nucl}_clean.pdb --all_histones histone-all.pdb --n_nucl ${n_nucl} --output_dir .
histone_dna_data_dir=$(pwd)
cd ..

# prepare rigid group file
cp ${lammps_data_main_dir}/single_nucleo/post_processing/groupsDefinition/group_rigid.txt .
python ${ca_sbm_3spn2c_openmm}/scripts/clean_group_rigid.py --input group_rigid.txt --n_nucl ${n_nucl}

# prepare dna sequence file
cp ${lammps_data_main_dir}/single_nucleo/buildDna/dnaSeq.txt .

# prepare smog files
mkdir -p smog-outputs
smog_output_dir=smog-outputs
cd smog-outputs
cp ${lammps_data_main_dir}/single_nucleo/sbm_scripts/single_nucl/smog2/smog.top .
cd ..

# build openmm system
python ${ca_sbm_3spn2c_openmm}/build_mono_fiber.py --env_main_dir ${ca_sbm_3spn2c_openmm} --n_nucl ${n_nucl} --histone_dna_data_dir ${histone_dna_data_dir} --group_rigid_txt_path group_rigid.txt --main_output_dir . --smog_output_dir ${smog_output_dir} --dna_seq_file dnaSeq.txt --temp 300 --salt 150 --platform CPU --mode default_exclude_CA_1_4 --periodic --cubic_box_length 1000

cd ..
