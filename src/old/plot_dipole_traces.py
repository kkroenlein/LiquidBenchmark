import chemistry
from gaff2xml import cirpy
import mdtraj as md
import pymbar
import os
import pandas as pd
import glob
from density_simulation_parameters import DATA_PATH

num_bootstrap = 100
fixed_block_length = 20  # 200 ps blocks for dielectric error bar block averaging.

prmtop_filenames = glob.glob(DATA_PATH + "/tleap/*.prmtop")
filename_munger = lambda filename: os.path.splitext(os.path.split(filename)[1])[0].split("_")

for prmtop_filename in prmtop_filenames[-5:]:
    cas, n_molecules, temperature = filename_munger(prmtop_filename)
    print(cas, temperature)
    dcd_filename = DATA_PATH + "/production/%s_%s_%s_production.dcd" % (cas, n_molecules, temperature)
    csv_filename = DATA_PATH + "/production/%s_%s_%s_production.csv" % (cas, n_molecules, temperature)
    try:
        traj = md.load(dcd_filename, top=prmtop_filename)
    except IOError:
        continue
    if traj.unitcell_lengths is None: continue
    rho = pd.read_csv(csv_filename)["Density (g/mL)"].values * 1000.  # g / mL -> kg /m3
    initial_traj_length = len(traj)
    initial_density_length = len(rho)
    [t0, g, Neff] = pymbar.timeseries.detectEquilibration(rho)
    mu = rho[t0:].mean()
    sigma = rho[t0:].std() * Neff ** -0.5
    prmtop = chemistry.load_file(prmtop_filename)
    charges = prmtop.to_dataframe().charge.values
    temperature = float(temperature)
    traj = traj[t0 * len(traj) / len(rho):]
    formula = cirpy.resolve(cas, "formula")
    dip = md.geometry.dipole_moments(traj, charges)    
    figure()
    plot(dip)
