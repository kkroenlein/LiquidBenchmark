import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("whitegrid")


expt = pd.read_csv("./tables/data_dielectric.csv")
expt["temperature"] = expt["Temperature, K"]


pred = pd.read_csv("./tables/predictions.csv")
pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates  # Should be fixed now, probably due to the CAS / name duplication issue found by Julie.
#expt = expt.groupby(["cas", "temperature"]).mean()  # Fix a couple of duplicates, not sure how they got there.
pred = pred.set_index(["cas", "temperature"])

pred["expt_density"] = expt["Mass density, kg/m3"]
pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]
pred["expt_density_std"] = expt["Mass density, kg/m3_std"]
pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_std"]


for (formula, grp) in pred.groupby("formula"):
    x, y = grp["density"], grp["expt_density"]
    xerr = grp["density_sigma"]
    yerr = grp["expt_density_std"].replace(np.nan, 0.0)
    plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label=formula)

plt.plot([600, 1400], [600, 1400], 'k')
plt.title("Density [kg / m^3]")
plt.xlim((600, 1400))
plt.ylim((600, 1400))
plt.xlabel("Predicted (GAFF)")
plt.ylabel("Experiment (ThermoML)")
plt.savefig("./manuscript/figures/densities_thermoml.pdf", bbox_inches=None)


yerr = pred["expt_dielectric_std"].replace(np.nan, 0.0)
xerr = pred["dielectric_sigma"].replace(np.nan, 0.0)

plt.figure()

plt.xlabel("Predicted (GAFF)")
plt.ylabel("Experiment (ThermoML)")
title("Static Dielectric Constant")

ticks = np.concatenate([np.arange(1, 10), 10 * np.arange(1, 10)])

xticks(ticks)
yticks(ticks)

plt.plot([1, 100], [1, 100], 'k')  # Guide
xscale('log')
yscale('log')


x, y = pred["dielectric"], pred["expt_dielectric"]
ols_model = sm.OLS(y, x)
ols_results = ols_model.fit()
r2 = ols_results.rsquared
#plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="GAFF (R^2 = %.3f)" % r2)
plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="GAFF")

xlim((1, 100))
ylim((1, 100))
plt.legend(loc=0)
plt.savefig("./manuscript/figures/dielectrics_thermoml_nocorr.pdf", bbox_inches=None)


x, y = pred["corrected_dielectric"], pred["expt_dielectric"]
ols_model = sm.OLS(y, x)
ols_results = ols_model.fit()
r2 = ols_results.rsquared
#plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="Corrected (R^2 = %.3f)" % r2)
plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="Corrected")

xlim((1, 100))
ylim((1, 100))
plt.legend(loc=0)
plt.savefig("./manuscript/figures/dielectrics_thermoml.pdf", bbox_inches=None)
