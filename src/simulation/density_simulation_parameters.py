"""
This file contains all simulation parameters used in this project.
"""
from simtk import unit as u

CUTOFF = 0.95 * u.nanometers

PRESSURE = 1.0 * u.atmospheres
BAROSTAT_FREQUENCY = 25

FRICTION = 1.0 / u.picoseconds
EQUIL_FRICTION = 5.0 / u.picoseconds


EQUIL_TIMESTEP = 1.0 * u.femtoseconds
TIMESTEP = 1.0 * u.femtoseconds


N_STEPS = 500000 # 0.5 ns (at a time)
N_EQUIL_STEPS = 5000000 # 5ns

OUTPUT_FREQUENCY = 10000 # 10ps
OUTPUT_DATA_FREQUENCY = 250 # 0.25ps

STD_ERROR_TOLERANCE = 0.0002 # g/mL


# DEBUG PARAMETERS
N_STEPS = 100000
N_EQUIL_STEPS = 50000
STD_ERROR_TOLERANCE = 0.004 # g/mL
