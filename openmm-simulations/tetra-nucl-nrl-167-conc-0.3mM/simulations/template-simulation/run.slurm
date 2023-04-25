#!/bin/bash
#SBATCH --job-name=md
#SBATCH --output=slurm-output/md-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=20GB
#SBATCH --time=2-00:00:00

module load cuda/10.0

snapshot=../../test-snapshots/snapshot0.pdb

echo "Run step1 relax"
python step1_relax.py --input_tetra_nucl_pdb ${snapshot} --output_dcd relax.dcd --output_interval 10000 --n_steps 200000

# for step2 and step3, input pdb is the structure including tetranucleosome and single nucleosomes, and this pdb file only provides the topology
# for step2, the initial coordinate is read from input_dcd
echo "Run step2 energy minimization"
python step2_min.py --input_pdb start.pdb --input_dcd relax.dcd --output_state min.xml

# for step3, the initial coordinate is read from input_state
echo "Run step3 production"
python step3_prod.py --input_pdb start.pdb --input_state min.xml --plumed plumed.txt --output_dcd prod.dcd --output_interval 10000 --n_steps 1000000

echo "job done"

