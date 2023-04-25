import numpy as np
import pandas as pd
import sys
import os
import argparse
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import warnings
warnings.filterwarnings("ignore")

sys.path.append('/home/gridsan/sliu/Projects/smog-3spn2-openmm')
from OpenSMOG3SPN2.forcefields.rigid import createRigidBodies
import OpenSMOG3SPN2.utils.helper_functions as helper_functions
from OpenSMOG3SPN2.utils.chromatin_helper_functions import get_chromatin_rigid_bodies
from OpenSMOG3SPN2.utils.insert import insert_molecules

parser = argparse.ArgumentParser()
parser.add_argument('--input_tetra_nucl_pdb', required=True, help='input tetranucleosome pdb path')
parser.add_argument('--output_dcd', default='relax.dcd', help='output dcd path')
parser.add_argument('--output_data', default='relax.log', help='output state data path')
parser.add_argument('--output_interval', type=int, default=100000, help='output interval')
parser.add_argument('--n_steps', type=int, default=1000000, help='simulation number of steps')
args = parser.parse_args()

build_dir = '../../build-system'
nonrigid_system_xml = f'{build_dir}/nonrigid_system.xml'

# prepare initial configuration
cg_single_nucl_pdb = f'{build_dir}/cg_single_nucl.pdb'
n_single_nucl = 26
box_a, box_b, box_c = 55.0, 55.0, 55.0
tetra_nucl_atoms = helper_functions.parse_pdb(args.input_tetra_nucl_pdb)
coord = tetra_nucl_atoms[['x', 'y', 'z']].to_numpy()
coord -= np.mean(coord, axis=0)
coord += 0.5*10*np.array([box_a, box_b, box_c])
tetra_nucl_atoms[['x', 'y', 'z']] = coord
helper_functions.write_pdb(tetra_nucl_atoms, 'aligned_tetra_nucl.pdb')
insert_molecules(cg_single_nucl_pdb, 'start.pdb', n_mol=n_single_nucl, existing_pdb='aligned_tetra_nucl.pdb', 
                 radius=1.0, box=[box_a, box_b, box_c])
top = app.PDBFile('start.pdb').getTopology()
n_atoms = top.getNumAtoms()
init_coord = app.PDBFile('start.pdb').getPositions()
rigid_coord = init_coord

# set rigid bodies
# fix the whole tetranucleosome, and for each single nucleosome, fix histone core with middle 73 bp core DNA
rigid_bodies = []
n_tetra_nucl_atoms = app.PDBFile(args.input_tetra_nucl_pdb).getTopology().getNumAtoms()
print(f'{n_tetra_nucl_atoms} atoms in tetranucleosome.')
rigid_bodies.append(list(range(n_tetra_nucl_atoms)))

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
simulation.context.setVelocitiesToTemperature(temperature)

dcd_reporter = app.DCDReporter(args.output_dcd, args.output_interval, enforcePeriodicBox=True)
state_reporter = app.StateDataReporter(args.output_data, args.output_interval, step=True, time=True, 
                                       potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, 
                                       speed=True)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)
simulation.step(args.n_steps)


