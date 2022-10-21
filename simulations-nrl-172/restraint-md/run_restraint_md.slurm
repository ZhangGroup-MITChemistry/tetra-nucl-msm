#!/bin/bash
#SBATCH --job-name=tetra
#SBATCH --array=2-200
#SBATCH --output=slurm-output/slurm-%A_%a.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=10GB
#SBATCH --time=12:00:00

module load cuda/10.0

job_id=${SLURM_ARRAY_TASK_ID}
curr=$(printf "%05d" ${job_id})

scripts_dir=/home/gridsan/sliu/Projects/tetra-nucl-msm/scripts
template_plumed_restraint=${scripts_dir}/template_plumed_restraint.txt

mkdir -p run${curr}
cd run${curr}
python ${scripts_dir}/prep_restraint_plumed.py --job_id ${job_id} --template_plumed_restraint ${template_plumed_restraint} --output_plumed_restraint plumed_restraint.txt

nrl=172
ca_sbm_3spn2c_openmm=/home/gridsan/sliu/Projects/chromatin-aggregation/CA_SBM_3SPN2C_OPENMM
system_xml=/home/gridsan/sliu/Projects/tetra-nucl-msm/openmm-structures/nrl-${nrl}-4mer/mode-default_exclude_CA_1_4-150mM-300K-PBC-box-1000.00nm-init-system-state/system.xml
state_xml=/home/gridsan/sliu/Projects/tetra-nucl-msm/openmm-structures/nrl-${nrl}-4mer/mode-default_exclude_CA_1_4-150mM-300K-PBC-box-1000.00nm-init-system-state/state.xml
cg_fiber_pdb=/home/gridsan/sliu/Projects/tetra-nucl-msm/openmm-structures/nrl-${nrl}-4mer/cg-fiber/cg_fiber.pdb
python ${ca_sbm_3spn2c_openmm}/run_system.py --env_main_dir ${ca_sbm_3spn2c_openmm} --system ${system_xml} --state ${state_xml} --pdb ${cg_fiber_pdb} --move_COM_to_box_center -o . --platform CUDA --precision mixed --temp 300 --simulation_mode const_temp --minimize --report_freq 20000 --steps 200000 --plumed plumed_restraint.txt


