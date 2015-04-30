"""
Builds a table of error by functional groups.  Not currently used in the paper because the number of unique chemicals is quite low.
"""
import sklearn.metrics, sklearn.cross_validation
import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess
import openmoltools

X = pd.read_csv("./tables/data_dielectric.csv")

data = []
for k, smiles in enumerate(X.smiles.unique()):
    filename = "./mol2/%d.mol2" % k
    if not os.path.exists(filename):
        oemol = openmoltools.openeye.smiles_to_oemol(smiles)
        oemol = openmoltools.openeye.generate_conformers(oemol, strictStereo=False)
        openmoltools.openeye.molecule_to_mol2(oemol, filename)
    checkmol = subprocess.check_output(["checkmol", filename])
    data.append(dict(smiles=smiles, checkmol=checkmol))

data = pd.DataFrame(data)
v = data.checkmol.map(lambda x: x.split("\n"))
groups = unique(sum(v))
groups = np.array([g for g in groups if g not in ["", "cation"]])


for group in groups:
    data[group] = data.checkmol.map(lambda x: group in x)


counts = data[groups].sum(0)
counts = pd.DataFrame(counts)
print(counts.to_latex())



expt = pd.read_csv("./tables/data_with_metadata.csv")
expt["temperature"] = expt["Temperature, K"]


pred = pd.read_csv("./tables/predictions.csv")
pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates  # Should be fixed now, probably due to the CAS / name duplication issue found by Julie.
#expt = expt.groupby(["cas", "temperature"]).mean()  # Fix a couple of duplicates, not sure how they got there.
pred = pred.set_index(["cas", "temperature"])

pred["expt_density"] = expt["Mass density, kg/m3"]
pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]
#pred["expt_density_std"] = expt["Mass density, kg/m3_std"]
pred["expt_density_std"] = expt["Mass density, kg/m3_uncertainty_bestguess"]
#pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_std"]
pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_uncertainty_bestguess"]



x, y = pred["density"], pred["expt_density"]
smiles_to_cas = pd.Series(X[["smiles", "cas"]].set_index("smiles").to_dict()["cas"])
cas_to_smiles = pd.Series(X[["smiles", "cas"]].set_index("cas").to_dict()["smiles"])
pred["smiles"] = pred.cas.apply(lambda x: cas_to_smiles[x])

pred["density_relative_error"] = (pred["density"] - pred["expt_density"]) / pred["expt_density"]
pred["density_relative_error_squared"] = pred.density_relative_error ** 2

D = data[["smiles"] + groups.tolist()].set_index("smiles")

for group in groups:
    pred[group] = pred.smiles.apply(lambda x: D[group][x])

baseline = pred.density_relative_error_squared.mean() ** 0.5
density_rms_pct_error = pd.DataFrame({group:pred.groupby(group).density_relative_error_squared.mean() ** 0.5 for group in groups}).T
density_rms_pct_error["ratio"] = density_rms_pct_error[True] / density_rms_pct_error[False]
density_rms_pct_error["over_baseline"] = density_rms_pct_error[True] / baseline
