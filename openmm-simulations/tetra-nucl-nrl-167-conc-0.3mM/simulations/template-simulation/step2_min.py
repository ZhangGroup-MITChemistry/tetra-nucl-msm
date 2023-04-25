import numpy as np
import pandas as pd
import sys
import os
import argparse
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj
import warnings
warnings.filterwarnings("ignore")

sys.path.append('/home/gridsan/sliu/Projects/smog-3spn2-openmm')
from OpenSMOG3SPN2.forcefields.rigid import createRigidBodies
from OpenSMOG3SPN2.utils.chromatin_helper_functions import get_chromatin_rigid_bodies

parser = argparse.ArgumentParser()
parser.add_argument('--input_pdb', required=True, help='input tetranucleosome pdb path')
parser.add_argument('--input_dcd', default='relax.dcd', help='input dcd path')
parser.add_argument('--output_state', default='min.xml', help='output state xml path')
args = parser.parse_args()

build_dir = '../../build-system'
nonrigid_system_xml = f'{build_dir}/nonrigid_system.xml'

cg_single_nucl_pdb = f'{build_dir}/cg_single_nucl.pdb'
n_single_nucl = 26
top = app.PDBFile(args.input_pdb).getTopology()
n_atoms = top.getNumAtoms()
traj = mdtraj.load_dcd(args.input_dcd, top=args.input_pdb)
init_coord = traj.xyz[-1]*unit.nanometer
rigid_coord = init_coord

# set rigid bodies
rigid_bodies = []
n_tetra_nucl_atoms = app.PDBFile(f'{build_dir}/cg_tetra_nucl.pdb').getTopology().getNumAtoms()
print(f'{n_tetra_nucl_atoms} atoms in tetranucleosome.')
tetra_nucl_rigid_bodies = get_chromatin_rigid_bodies(n_nucl=4, nrl=172)
rigid_bodies += tetra_nucl_rigid_bodies

single_nucl_rigid_body = np.array(get_chromatin_rigid_bodies(n_nucl=1, nrl=147)[0])
n_atoms_per_single_nucl = app.PDBFile(cg_single_nucl_pdb).getTopology().getNumAtoms()
print(f'{n_atoms_per_single_nucl} atoms in each single nucleosome')
assert n_atoms == n_tetra_nucl_atoms + n_single_nucl*n_atoms_per_single_nucl
for i in range(n_single_nucl):
    rigid_bodies.append((single_nucl_rigid_body + n_tetra_nucl_atoms + i*n_atoms_per_single_nucl).tolist())

with open(nonrigid_system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
createRigidBodies(system, rigid_coord, rigid_bodies)

temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
platform_name = 'CUDA'
platform = mm.Platform.getPlatformByName(platform_name)
properties = {'Precision': 'mixed'}
simulation = app.Simulation(top, system, integrator, platform, properties)
simulation.context.setPositions(init_coord)
simulation.minimizeEnergy()
state = simulation.context.getState(getPositions=True)
with open(args.output_state, 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))


