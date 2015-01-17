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




# HACK TO PLOT DATA THAT LACKS ERRORBARS!
pred.expt_density_std = pred.expt_density_std.fillna(0.0)

g = sns.FacetGrid(pred, col="name", col_wrap=6, size=3.5)
g.map(plt.errorbar, "temperature", "expt_density", "expt_density_std", fmt='.', color='b', label="Expt")
g.map(plt.errorbar, "temperature", "density", "density_sigma", fmt='.', color='g', label="MD")
g.set_ylabels("Density [kg / m^3")
legend(loc=4)
plt.savefig("./manuscript/figures/densities_versus_temperature_all.pdf", bbox_inches=None)


# HACK TO PLOT DATA THAT LACKS ERRORBARS!
pred.expt_dielectric_std = pred.expt_dielectric_std.fillna(0.0)
# NEED TO DO ERROR ESTIMATES ON EPSILON
pred["corrected_dielectric_sigma"] = 0.0
pred["dielectric_sigma"] = 0.0
# End HACKERY

g = sns.FacetGrid(pred, col="name", col_wrap=6, size=3.5)
g.set(yscale="log")
g.map(plt.errorbar, "temperature", "expt_dielectric", "expt_dielectric_std", fmt='.', color='b', label="Expt")
g.map(plt.errorbar, "temperature", "dielectric", "dielectric_sigma", fmt='.', color='g', label="MD")
g.map(plt.errorbar, "temperature", "corrected_dielectric", "corrected_dielectric_sigma", fmt='.', color='r', label="Corr.")
g.set(yticks=[0.5, 0,6, 0.7, 0.8, 0.9] + range(1, 10) + range(10, 100, 10) + range(100, 1000, 100))
g.set(ylim=(0.5, 200))
g.set_ylabels("Dielectric")
legend(loc=4)
plt.savefig("./manuscript/figures/dielectric_versus_temperature_all.pdf", bbox_inches=None)
