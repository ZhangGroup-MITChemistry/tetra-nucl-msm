#!/bin/bash
#SBATCH --job-name=tetra
#SBATCH --array=0-1000%10
#SBATCH --output=slurm-output/slurm-%A_%a.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=10GB
#SBATCH --time=12:00:00

module load cuda/10.0

# set job id
job_id=${SLURM_ARRAY_TASK_ID}
curr=$(printf "%05d" ${job_id})

python ${ca_sbm_main_dir}/run_system.py --system ${system_xml} --pdb ${pdb} --state ${state_xml} --coord input_coord.txt --move_COM_to_box_center --output_dir . --platform CUDA --precision mixed --minimize --steps ${n_steps} --report_freq ${freq} --plumed ${plumed} --temp ${temp} --integrator NoseHoover --collision 1.0 > simulation_output.txt


echo "job done"

