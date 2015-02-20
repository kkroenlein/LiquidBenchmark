import itertools
import pandas as pd
import mdtraj as md
import glob
import simtk.unit as u
import numpy as np

"""
Calculate the radius and dipole moment of each molecule in FreeSolv.
Then calculate the predicted error in solvation free energy in cyclohexane
that results from the missing polarization physics.

"""

ke = 8.98755E9 * u.meter / u.farad  # Coulomb constant (http://en.wikipedia.org/wiki/Coulomb%27s_constant)
epsilon = 2.023  # For cyclohexane: http://en.wikipedia.org/wiki/Cyclohexane_%28data_page%29
freesolv_path = "/home/kyleb/src/choderalab/FreeSolv/mol2files_gaff/"
filenames = glob.glob(freesolv_path + "*.mol2")

data = []
for filename in filenames:
    t = md.load(filename)
    atom_pairs = np.array(list(itertools.product(np.arange(t.n_atoms), np.arange(t.n_atoms))))
    r = 0.5 * md.compute_distances(t, atom_pairs).max() * u.nanometers
    frame, bonds = md.formats.mol2.mol2_to_dataframes(filename)
    q = frame["charge"].values
    mu = md.geometry.dipole_moments(t, q)
    mu = np.linalg.norm(mu) * u.nanometers * u.elementary_charge
    current_data = dict(filename=filename, mu=mu, r=r)
    data.append(current_data)

data = pd.DataFrame(data)
data["onsager_prefactor"] = -(data.mu ** 2) / (data.r ** 3) * ke

# Need to fix some weirdness going beween coulomb, elementary charge, and farad in the openmm units package.
unit_fixer = lambda x: (x / (u.coulomb ** 2 / u.farad) * u.joule * u.AVOGADRO_CONSTANT_NA) / u.kilocalories_per_mole
data["onsager_prefactor"] = data["onsager_prefactor"].apply(unit_fixer)

eps_term = lambda x: (x - 1.) / (2. * x + 1)
data["delta_G"] = (eps_term(epsilon) - eps_term(1.0)) * data["onsager_prefactor"]

data.delta_G.mean(), data.delta_G.std()
