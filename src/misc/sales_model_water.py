import polarizability
import simtk.unit as u
import numpy as np

alpha = (0.32 + 0.17 * 2 + 0.57) * u.angstrom ** 3
alpha1 = polarizability.polarizability_from_formula("H2O")
alpha, alpha1, alpha - alpha1

delta = 4 * pi * 55.4 * u.molar * alpha
delta = delta * u.AVOGADRO_CONSTANT_NA

r = polarizability.dielectric_correction_from_formula("H2O", 1000. * (u.gram / u.milliliter))
r, delta, r - delta

