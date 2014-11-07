import statsmodels.formula.api as sm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import simtk.unit as u
import polarizability

metadata = pd.read_table("./virtualchemistry.tab").set_index("cas")

data = pd.read_csv("./vchem.csv").set_index("cas")
data["formula"] = metadata.formula
data["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in data[["formula", "density"]].iterrows()))
data["corrected"] = data.expt - data.polcorr
data["gaffdiff"] = data.expt - data.gaff
#result = sm.ols(formula="gaffdiff ~ polcorr", data=data).fit()
data["optcorr"] = -1.369696 + 4.626496 * data.polcorr
data["optcorrected"] = data.expt - data.optcorr



g = sns.lmplot("gaff", "expt", data, fit_reg=False)
plt.plot([1, 100], [1, 100], 'k')  # Guide
g.set(xscale='log', xlim=(1, 100), ylim=(1, 100))
g.set(yscale='log')


g = sns.lmplot("gaff", "corrected", data, fit_reg=False)
plt.plot([1, 100], [1, 100], 'k')  # Guide
g.set(xscale='log', xlim=(1, 100), ylim=(1, 100))
g.set(yscale='log')

g = sns.lmplot("gaff", "optcorrected", data, fit_reg=False)
plt.plot([1, 100], [1, 100], 'k')  # Guide
g.set(xscale='log', xlim=(1, 100), ylim=(1, 100))
g.set(yscale='log')

figure()
g = sns.residplot("gaff", "expt", data)
g.set(xscale='log', xlim=(1, 100))


figure()
g = sns.residplot("gaff", "corrected", data)
g.set(xscale='log', xlim=(1, 100))


figure()
g = sns.residplot("gaff", "optcorrected", data)
g.set(xscale='log', xlim=(1, 100))




result = sm.ols(formula="gaffdiff ~ polcorr", data=data).fit()
