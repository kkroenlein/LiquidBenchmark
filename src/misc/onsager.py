import simtk.unit as u
import numpy as np

debye = 3.33564E-30 * u.coulomb * u.meter  # simtk.unit debye is messed up, using http://en.wikipedia.org/wiki/Debye
ke = 8.98755E9 * u.meter / u.farad  # Coulomb constant (http://en.wikipedia.org/wiki/Coulomb%27s_constant)

# Function to calculate microscopic transfer free energy due to electrostatics in onsager model
f0 = lambda epsilon, mu, d: ke * mu ** 2 * d ** -3 * (epsilon - 1) / (2. * epsilon + 1.)

# Function to convert to kcal / mol
f = lambda epsilon, mu, d: (f0(epsilon, mu, d) / (u.coulomb ** 2 / u.farad) * u.joule * u.AVOGADRO_CONSTANT_NA) / u.kilocalories_per_mole

mu = 2.0 * debye  # Similar to water.  http://en.wikipedia.org/wiki/Properties_of_water
d = 2.0 * u.angstrom  # Made up, but similar to water
f(2.2, mu, d) - f(1.0, mu, d)  # CCl4
# 


mu = 2.2 * debye  # Similar to water.  http://en.wikipedia.org/wiki/Properties_of_water
d = 1.93 * u.angstrom  # Made up, but similar to water
f(2.2, mu, d) - f(1.0, mu, d)  # CCl4
# 
