from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import openmmtools

cas = "126492-54-4"

#integrator_type = "langevin"
#integrator_type = "langevin1fs"
#integrator_type = "ghmc"
integrator_type = "hmc"

steps_per_hmc = 10

#output_frequency = 250
output_frequency = 500 if integrator_type != "hmc" else 500 / steps_per_hmc
dcd_frequency = output_frequency * 20

timestep = 2.0 * u.femtoseconds
friction = 1.0 / u.picoseconds
temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres
barostat_frequency = 25
n_steps = 500000000
cutoff = 1.0 * u.nanometers

ffxml_filename = "%s.xml" % cas
ff = app.ForceField(ffxml_filename)

dcd_filename = "./production/production_%s.dcd" % integrator_type
log_filename = "./production/production_%s.log" % integrator_type

pdb = app.PDBFile("./equil/equil_final_step.pdb")

topology = pdb.topology
positions = pdb.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrators = {"langevin":mm.LangevinIntegrator(temperature, friction, timestep), 
"ghmc":openmmtools.integrators.GHMCIntegrator(temperature, friction, timestep),
"hmc":openmmtools.integrators.HMCIntegrator(temperature, steps_per_hmc, timestep),
"langevin1fs":mm.LangevinIntegrator(temperature, friction, timestep / 2.),
}

integrator = integrators[integrator_type]

system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

simulation = app.Simulation(topology, system, integrator)

simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(app.DCDReporter(dcd_filename, dcd_frequency))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, density=True))
print(integrator_type, output_frequency)
simulation.step(n_steps)
