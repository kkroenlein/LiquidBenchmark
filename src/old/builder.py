import numpy as np
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

def build_simulation(traj, forcefield):

    topology = traj.top.to_openmm(traj[0])
    positions = traj.openmm_positions(0)

    system = forcefield.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9*unit.nanometers)

    charges = []
    for force in system.getForces():
        if hasattr(force, "getParticleParameters"):
            for k in range(force.getNumParticles()):
                q, s, e = force.getParticleParameters(k)
                charges.append(q / unit.elementary_charge)
    
    charges = np.array(charges)
    
    return system, charges

    integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds)
    system.addForce(mm.MonteCarloBarostat(1 * unit.atmospheres, 300 * unit.kelvin, 25))

    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    return simulation
