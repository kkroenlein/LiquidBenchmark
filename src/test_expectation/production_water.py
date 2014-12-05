from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import openmmtools

cas = "tip3p"

#integrator_type = "langevin"
integrator_type = "langevin1fs"
#integrator_type = "ghmc"
#integrator_type = "hmc"

output_frequency = 250
dcd_frequency = 5000

steps_per_hmc = 250

timestep = 2.0 * u.femtoseconds
friction = 1.0 / u.picoseconds
temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres
barostat_frequency = 25
n_steps = 500000000
cutoff = 1.0 * u.nanometers
output_frequency = 500

ffxml_filename = "%s.xml" % cas
ff = app.ForceField(ffxml_filename)

dcd_filename = "./water/production_%s.dcd" % integrator_type
log_filename = "./water/production_%s.log" % integrator_type

pdb = app.PDBFile("../tip3p.pdb")

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
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(app.DCDReporter(dcd_filename, dcd_frequency))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, density=True))
print(integrator_type)
simulation.step(n_steps)
