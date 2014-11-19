import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("whitegrid")


expt = pd.read_csv("./data_dielectric.csv")
expt["temperature"] = expt["Temperature, K"]


pred = pd.read_csv("./predictions.csv")
pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

#expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates
expt = expt.groupby(["cas", "temperature"]).mean()  # Fix a couple of duplicates, not sure how they got there.
pred = pred.set_index(["cas", "temperature"])

pred["expt_density"] = expt["Mass density, kg/m3"]
pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]



x, y = pred["density"], pred["expt_density"]
plt.plot(x, y, 'o')

plt.title("Density [kg / m^3]")
plt.xlabel("Predicted")
plt.ylabel("Experiment")
plt.savefig("./manuscript/figures/densities_thermoml.pdf", bbox_inches=None)


plt.figure()
x, y = pred["dielectric"], pred["expt_dielectric"]
ols_model = sm.OLS(y, x)
ols_results = ols_model.fit()
r2 = ols_results.rsquared
plot(x, y, 'o', label="GAFF (R^2 = %.3f)" % r2)

x, y = pred["corrected_dielectric"], pred["expt_dielectric"]
ols_model = sm.OLS(y, x)
ols_results = ols_model.fit()
r2 = ols_results.rsquared
plot(x, y, 'o', label="Corrected (R^2 = %.3f)" % r2)

plt.plot([1, 100], [1, 100], 'k')  # Guide
xscale('log')
yscale('log')
xlim((1, 100))
ylim((1, 100))
plt.legend(loc=0)

ticks = np.concatenate([np.arange(1, 10), 10 * np.arange(1, 10)])

xticks(ticks)
yticks(ticks)

xlabel("Predicted")
ylabel("Experiment")
title("Static Dielectric Constant")
plt.savefig("./manuscript/figures/dielectrics_thermoml.pdf", bbox_inches=None)
