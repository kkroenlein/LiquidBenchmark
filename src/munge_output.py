from thermopyl import cirpy
from simtk.openmm import app
import builder
import mdtraj as md
import pymbar
import scipy.interpolate
import os
import pandas as pd
import glob

filenames = glob.glob("./data/equil/*.pdb")
filename_munger = lambda filename: os.path.splitext(os.path.split(filename)[1])[0].split("_")
data = []
for pdb_filename in filenames:
    cas, n_molecules, temperature, stage = filename_munger(pdb_filename)
    print(cas, temperature)
    dcd_filename = "./data/production/%s_%s_%s_production.dcd" % (cas, n_molecules, temperature)
    csv_filename = "./data/production/%s_%s_%s_production.csv" % (cas, n_molecules, temperature)
    try:
        traj = md.load(dcd_filename, top=pdb_filename)
    except IOError:
        continue
    if traj.unitcell_lengths is None: continue
    rho = pd.read_csv(csv_filename)["Density (g/mL)"].values
    #rho = md.geometry.density(traj)
    [t0, g, Neff] = pymbar.timeseries.detectEquilibration(rho)
    mu = rho[t0:].mean()
    sigma = rho[t0:].std() * Neff ** -0.5
    forcefield = app.ForceField("./data/ffxml/%s.xml" % cas)
    system, charges = builder.build_simulation(traj, forcefield)
    temperature = float(temperature)
    dielectric = md.geometry.static_dielectric(traj, charges, temperature)
    formula = cirpy.resolve(cas, "formula")
    data.append(dict(cas=cas, temperature=temperature, density=mu, density_sigma=sigma, Neff=Neff, dielectric=dielectric, formula=formula))
    print(data[-1])

data = pd.DataFrame(data)

data.to_csv("./tables/predictions.csv")
