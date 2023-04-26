#!/bin/bash
#SBATCH --job-name=md
#SBATCH --output=slurm-output/md-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=20GB
#SBATCH --time=2-00:00:00

module load cuda/10.0

# energy.csv and backup_system_energy.csv should be almost identical

python rerun.py ../build-system/nonrigid_system.xml energy.csv

python rerun.py ../build-system/backup/nonrigid_system.xml backup_system_energy.csv

echo "job done"

