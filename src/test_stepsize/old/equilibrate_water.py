import openmoltools
import mdtraj as md
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import openmmtools

temperature = 300. * u.kelvin
equil_friction = 5 / u.picoseconds
equil_timestep = 1 * u.femtoseconds
pressure = 1.0 * u.atmospheres
barostat_frequency = 25
discard_steps = 10000
n_steps = 500000
cutoff = 1.0 * u.nanometers
output_frequency = 500

n_monomers = 500
cas = "tip3p.xml"

ffxml_filename = "%s.xml" % cas
ff = app.ForceField(ffxml_filename)

pdb_filename = "./%s.pdb" % cas
box_pdb_filename = "./box.pdb"

monomer_pdb_filenames = [pdb_filename]
packed_trj = openmoltools.packmol.pack_box(monomer_pdb_filenames, [n_monomers])
packed_trj.save(box_pdb_filename)

out_pdb_filename = "./equil/equil.pdb"
final_step_pdb_filename = "./equil/equil_final_step.pdb"
dcd_filename = "./equil/equil.dcd"
log_filename = "./equil/equil.log"

topology = packed_trj.top.to_openmm(packed_trj)
positions = packed_trj.openmm_positions(0)

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, equil_friction, equil_timestep / 10.)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating...')

simulation.step(discard_steps)  # Don't even save the first XXX ps
integrator.setStepSize(equil_timestep)

simulation.reporters.append(app.DCDReporter(dcd_filename, n_steps - 1))
simulation.reporters.append(app.PDBReporter(out_pdb_filename, n_steps - 1))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, density=True))
simulation.step(n_steps)

del simulation
del system
t = md.load(dcd_filename, top=out_pdb_filename)
t0 = t[-1]
t0.unitcell_lengths = t.unitcell_lengths.mean(0)
t0.save(out_pdb_filename)

del t
del t0
t = md.load(dcd_filename, top=out_pdb_filename)[-1]
t.save(final_step_pdb_filename)
