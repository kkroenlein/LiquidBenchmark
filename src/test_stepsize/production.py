from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

cas = "126492-54-4"

timestep_factor = 1.

base_filename = "langevin%0.3fs" % 2.0 / timestep_factor

output_frequency = 500 * timestep_factor
dcd_frequency = output_frequency * 20

friction = 1.0 / u.picoseconds
temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres
barostat_frequency = 25
n_steps = 500000000
cutoff = 1.0 * u.nanometers

ffxml_filename = "%s.xml" % cas
ff = app.ForceField(ffxml_filename)

dcd_filename = "./production/production_%s.dcd" % base_filename
log_filename = "./production/production_%s.log" % base_filename

pdb = app.PDBFile("./equil/equil_final_step.pdb")

topology = pdb.topology
positions = pdb.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator mm.LangevinIntegrator(temperature, friction, 2.0 * u.femtoseconds / timestep_factor)

system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

simulation = app.Simulation(topology, system, integrator)

simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(app.DCDReporter(dcd_filename, dcd_frequency))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, density=True))
print(integrator_type, output_frequency)
simulation.step(n_steps)
