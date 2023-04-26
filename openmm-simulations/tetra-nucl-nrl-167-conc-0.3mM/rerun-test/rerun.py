import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj
import warnings
warnings.filterwarnings("ignore")

nonrigid_system_xml = sys.argv[1]
output_csv = sys.argv[2]

with open(nonrigid_system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())

top = app.PDBFile('start.pdb').getTopology()
#init_coord = app.PDBFile('start.pdb').getPositions()
temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
platform_name = 'CUDA'
platform = mm.Platform.getPlatformByName(platform_name)
properties = {'Precision': 'mixed'}
simulation = app.Simulation(top, system, integrator, platform, properties)

traj = mdtraj.load_dcd('traj.dcd', 'start.pdb')
df = pd.DataFrame(columns=[f'group {i}' for i in range(1, 13)])
for i in range(traj.xyz.shape[0]):
    simulation.context.setPositions(traj.xyz[i])
    row = []
    for j in range(1, 13):
        state = simulation.context.getState(getEnergy=True, groups={j})
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        row.append(energy)
    df.loc[len(df.index)] = row

#print(df)
df.round(6).to_csv(output_csv, index=False)

