from rdkit import Chem
import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("whitegrid")

d = pd.read_table("/home/kyleb/src/choderalab/oeante/oeante/gaff/gaffsmarts.txt", names=["smarts", "atype"], sep=r"\s*", comment='"""', engine="python", index_col=0)

expt = pd.read_csv("./data_dielectric.csv")
expt["temperature"] = expt["Temperature, K"]

m = Chem.MolFromSmiles('c1ccccc1O')

results = {}
for smiles in expt.smiles.unique():
    m = Chem.MolFromSmiles(smiles)
    results[smiles] = {}
    for (smarts, gaff) in d.itertuples():
        patt = Chem.MolFromSmarts(smarts)
        has_match = m.HasSubstructMatch(patt)
        matches = m.GetSubstructMatches(patt)
        results[smiles][gaff] = len(matches)

results = pd.DataFrame(results).T
results.index.name = "smiles"
results[results == 0] = np.nan
results.dropna(axis=1, how='all', inplace=True)
results.replace(np.nan, 0, inplace=True)
results[results > 0] = 1  # Binarize

pred = pd.read_csv("./predictions.csv")
pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

#expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates
smiles = expt.groupby(["cas", "temperature"]).first().smiles
expt = expt.groupby(["cas", "temperature"]).mean()  # Fix a couple of duplicates, not sure how they got there.

pred = pred.set_index(["cas", "temperature"])
pred["smiles"] = smiles

pred["expt_density"] = expt["Mass density, kg/m3"]
pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]

pred["logerr_density"] = np.log(pred.expt_density / pred.density)
pred["logerr_dielectric"] = np.log(pred.expt_dielectric / pred.corrected_dielectric)
pred["logerr_dielectric_raw"] = np.log(pred.expt_dielectric / pred.dielectric)

full = pd.merge(pred, results.reset_index(), on="smiles")



for measurement in ["logerr_density", "logerr_dielectric", "logerr_dielectric_raw"]:
    plt.figure()
    data = pd.DataFrame(dict((key, full.groupby(key)[measurement].mean()) for key in results.columns)).T
    data.index.name = "atom"
    data.columns = [False, True]
    data.T.boxplot()
    title(measurement)
    plt.savefig("/home/kyleb/src/kyleabeauchamp/LiquidBenchmark/manuscript/figures/functional_group_%s.pdf" % measurement, bbox_inches=None)
