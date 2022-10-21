#!/bin/bash
#SBATCH --job-name=submit
#SBATCH --output=slurm-output/submit_restraint.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=sched_mit_binz,newnodes,sched_mit_hill
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --exclude=node[445,455,457,462]
#SBATCH --export=ALL
#SBATCH --constraint=centos7

module purge
module load gcc
module load engaging/openmpi/1.8.8

scripts_dir=/nfs/pool002/users/smliu/tetra-nucl-msm/scripts
template_plumed_restraint=${scripts_dir}/template_plumed_restraint.txt

start_job_id=2
end_job_id=499

for (( i=${start_job_id}; i<=${end_job_id}; i++ ))
do 
    curr=$(printf "%05d" ${i})
    if [ ! -d restraint-md/run${curr} ]
    then
        cp -r restraint-md/template-run restraint-md/run${curr}
        cd restraint-md/run${curr}
        python ${scripts_dir}/prep_restraint_plumed.py --job_id ${i} --template_plumed_restraint ${template_plumed_restraint} --output_plumed_restraint plumed.txt
        source cmd.pre_submitJob.sh
        cd ../..
    fi
done

