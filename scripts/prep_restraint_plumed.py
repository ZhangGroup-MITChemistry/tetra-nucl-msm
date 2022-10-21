import numpy as np
import sys
import os
import argparse

# use this script to prepare plumed.txt file for restraint md
# we need to insert the restraint bias center

parser = argparse.ArgumentParser()
parser.add_argument('--job_id', type=int, required=True, help='input job id')
parser.add_argument('--template_plumed_restraint', default='template_plumed_restraint.txt', help='template plumed file')
parser.add_argument('--output_plumed_restraint', default='plumed_restraint.txt', help='output plumed file')
args = parser.parse_args()

dist_txt = '/nfs/pool002/users/smliu/tetra-nucl-msm/prep-init-snapshots/selected_sample_dist.txt'
dist_data = np.loadtxt(dist_txt, skiprows=1)

job_id = args.job_id
if job_id <= dist_data.shape[0] - 1:
    d = dist_data[job_id]
    d = d.tolist()
    d = [str(each) for each in d]

with open(args.template_plumed_restraint, 'r') as input_reader:
    template_plumed_restraint = input_reader.read()

plumed_restraint = template_plumed_restraint.replace('INPUT_AT', ','.join(d))

with open(args.output_plumed_restraint, 'w') as output_writer:
    output_writer.write(plumed_restraint)


