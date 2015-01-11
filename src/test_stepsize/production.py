from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

cas = "126492-54-4"

timestep_factor = 1.

base_filename = "langevin%0.1fs" % (2.0 / timestep_factor)
timestep = 2.0 * u.femtoseconds / timestep_factor

output_frequency = int(500 * timestep_factor)
dcd_frequency = int(output_frequency * 20)
barostat_frequency = int(25 * timestep_factor)
n_steps = int(500000000 * timestep_factor)

friction = 1.0 / u.picoseconds
temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres

cutoff = 1.0 * u.nanometers

ffxml_filename = "%s.xml" % cas
ff = app.ForceField(ffxml_filename)

dcd_filename = "./production/production_%s.dcd" % base_filename
log_filename = "./production/production_%s.log" % base_filename

pdb = app.PDBFile("./equil/equil_final_step.pdb")

topology = pdb.topology
positions = pdb.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, friction, timestep)

system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

simulation = app.Simulation(topology, system, integrator)

simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(app.DCDReporter(dcd_filename, dcd_frequency))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, density=True))
simulation.step(n_steps)
