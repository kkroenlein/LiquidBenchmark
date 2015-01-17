import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("whitegrid")


expt = pd.read_csv("./tables/data_with_metadata.csv")
expt["temperature"] = expt["Temperature, K"]


pred = pd.read_csv("./tables/predictions.csv")
pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates  # Should be fixed now, probably due to the CAS / name duplication issue found by Julie.
pred = pred.set_index(["cas", "temperature"])

pred["name"] = expt["components"]
pred["expt_density"] = expt["Mass density, kg/m3"]
pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]
pred["expt_density_std"] = expt["Mass density, kg/m3_uncertainty_bestguess"]
pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_uncertainty_bestguess"]

pred = pred.reset_index()


MIN_LENGTH = 0 
for (cas, grp) in pred.groupby("cas"):
    name = grp.name.iloc[0]
    temperature = grp["temperature"].values
    yhat, ytrue = grp["density"].values, grp["expt_density"].values
    yhat_err = grp["density_sigma"].values
    ytrue_err = grp["expt_density_std"].replace(np.nan, 0.0).values    
    if len(grp) > MIN_LENGTH:
        plt.figure()
    else:
        plt.figure(0)
    plt.errorbar(temperature, yhat, yerr=yhat_err, fmt='.', label="MD")
    plt.errorbar(temperature, ytrue, yerr=ytrue_err, fmt='.', label="Expt.")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Density [kg / m^3]")
    plt.legend(loc=0)
    if len(grp) > MIN_LENGTH:
        plt.title("%s (%s)" % (name, cas))
        #plt.savefig("./manuscript/figures/densities_versus_temperature_%s.pdf" % name, bbox_inches=None)




yerr = pred["expt_dielectric_std"].replace(np.nan, 0.0)
plt.figure()
x, y = pred["dielectric"], pred["expt_dielectric"]
ols_model = sm.OLS(y, x)
ols_results = ols_model.fit()
r2 = ols_results.rsquared
plt.errorbar(x, y, yerr=yerr, fmt='.', label="GAFF (R^2 = %.3f)" % r2)

x, y = pred["corrected_dielectric"], pred["expt_dielectric"]
ols_model = sm.OLS(y, x)
ols_results = ols_model.fit()
r2 = ols_results.rsquared
plt.errorbar(x, y, yerr=yerr, fmt='.', label="Corrected (R^2 = %.3f)" % r2)

plt.plot([1, 100], [1, 100], 'k')  # Guide
xscale('log')
yscale('log')
xlim((1, 100))
ylim((1, 100))
plt.legend(loc=0)

ticks = np.concatenate([np.arange(1, 10), 10 * np.arange(1, 10)])

xticks(ticks)
yticks(ticks)

plt.xlabel("Predicted (GAFF)")
plt.ylabel("Experiment (ThermoML)")
title("Static Dielectric Constant")
plt.savefig("./manuscript/figures/dielectrics_thermoml.pdf", bbox_inches=None)
