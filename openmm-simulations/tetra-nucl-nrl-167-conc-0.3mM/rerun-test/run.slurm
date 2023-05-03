#!/bin/bash
#SBATCH --job-name=rerun
#SBATCH --output=slurm-output/rerun-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=20GB
#SBATCH --time=2-00:00:00

module load cuda/10.0

cp ../build-system/start.pdb .

python rerun.py ../build-system/nonrigid_system.xml energy.csv

python rerun.py ../build-system/backup/nonrigid_system.xml backup_system_energy.csv

echo "job done"

