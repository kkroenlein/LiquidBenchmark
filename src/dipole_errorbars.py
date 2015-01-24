import mdtraj as md
import numpy as np
import mdtraj.utils.unit.unit_definitions as u

kB = 1.3806488E-23 * u.joule / u.kelvin
epsilon0 = 8.854187817E-12 * u.farad / u.meter
gas_constant = 8.3144621 * u.joule / u.kelvin / u.mole


def dipole_moment_errorbars():
    """Modified from mdtraj.geometry.static_dielectic()."""
    moments = md.geometry.dipole_moments(traj, charges)
    
    mu = moments.mean(0)  # Mean over frames
    
    subtracted = moments - mu
    
    dipole_variance = (subtracted * subtracted).sum(-1).mean(0) * (u.elementary_charge * u.nanometers) ** 2.  # <M*M> - <M>*<M> = <(M - <M>) * (M - <M>)>

    volume = traj.unitcell_volumes.mean() * u.nanometers ** 3.  # Average box volume of trajectory
    
    static_dielectric_sigma = dipole_variance / (3 * kB * temperature * volume * epsilon0)  # Eq. 7 of Derivation of an improved simple point charge model for liquid water: SPC/A and SPC/L 



def bootstrap(traj, charges, temperature, block_length):
    n = traj.n_frames / block_length
    indices = np.array_split(np.arange(traj.n_frames), n)
    epsilon = np.zeros(n)
    for k, ind in enumerate(indices):
        t = traj[ind]
        epsilon[k] = md.geometry.static_dielectric(t, charges, temperature)
    return epsilon, epsilon.std() * n ** -0.5
