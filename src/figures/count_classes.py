import re
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import glob
from thermopyl import thermoml_lib, cirpy

data = pd.read_hdf("./data.h5", 'data')

# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!
bad_filenames = ["./10.1016/j.fluid.2013.12.014.xml"]
data = data[~data.filename.isin(bad_filenames)]
# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!

experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]  # , "Isothermal compressibility, 1/kPa", "Isobaric coefficient of expansion, 1/K"]

ind_list = [data[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
X = data.ix[ind]

name_to_formula = pd.read_hdf("./compound_name_to_formula.h5", 'data')
name_to_formula = name_to_formula.dropna()

X["n_components"] = X.components.apply(lambda x: len(x.split("__")))
X = X[X.n_components == 1]
X.dropna(axis=1, how='all', inplace=True)

counts_data = {}
counts_data["0.  Single Component"] = X.count()[experiments]

X["formula"] = X.components.apply(lambda chemical: name_to_formula[chemical])

heavy_atoms = ["N", "C", "O", "S", "Cl", "Br", "F"]
desired_atoms = ["H"] + heavy_atoms

X["n_atoms"] = X.formula.apply(lambda formula_string : thermoml_lib.count_atoms(formula_string))
X["n_heavy_atoms"] = X.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, heavy_atoms))
X["n_desired_atoms"] = X.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, desired_atoms))
X["n_other_atoms"] = X.n_atoms - X.n_desired_atoms

X = X[X.n_other_atoms == 0]

counts_data["1.  Druglike Elements"] = X.count()[experiments]

X = X[X.n_heavy_atoms > 0]
X = X[X.n_heavy_atoms <= 10]
X.dropna(axis=1, how='all', inplace=True)

counts_data["2.  Heavy Atoms"] = X.count()[experiments]

X["smiles"] = X.components.apply(lambda x: cirpy.resolve(x, "smiles"))  # This should be cached via sklearn.
X = X[X.smiles != None]
X = X.ix[X.smiles.dropna().index]
    
X["cas"] = X.components.apply(lambda x: thermoml_lib.get_first_entry(cirpy.resolve(x, "cas")))  # This should be cached via sklearn.
X = X[X.cas != None]
X = X.ix[X.cas.dropna().index]

X = X[X["Temperature, K"] > 270]
X = X[X["Temperature, K"] < 330]

counts_data["3.  Temperature"] = X.count()[experiments]

X = X[X["Pressure, kPa"] > 100.]
X = X[X["Pressure, kPa"] < 102.]

counts_data["4.  Pressure"] = X.count()[experiments]

X.dropna(axis=1, how='all', inplace=True)

X["Pressure, kPa"] = 101.325  # Assume everything within range is comparable.  
X["Temperature, K"] = X["Temperature, K"].apply(lambda x: x.round(1))  # Round at the 0.1 digit.  


mu = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].mean()
sigma = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].std().dropna()

counts_data["5.  Aggregate T, P"] = mu.count()[experiments]

counts_data = pd.DataFrame(counts_data).T

q = mu.reset_index()
q = q.ix[q[experiments].dropna().index]
counts_data.ix["6.  Density+Dielectric"] = len(q)

print counts_data.to_latex()

plt.figure()
X["Temperature, K"].hist()
plt.title("ThermoML Density Data")
plt.ylabel("Number of Measurements")
plt.xlabel("Temperature [K]")

plt.savefig("/home/kyleb/src/kyleabeauchamp/LiquidBenchmark/manuscript/figures/thermoml_density_histogram.pdf", bbox_inches=None)
