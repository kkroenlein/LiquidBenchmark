import matplotlib.pyplot as plt
import pandas as pd

experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]

data = pd.read_csv("./tables/data_with_metadata.csv")

x = data["Mass density, kg/m3_uncertainty_author"]
y = data["Mass density, kg/m3_uncertainty_std"]
ind = x.dropna().index.intersection(y.dropna().index)
x, y = x[ind], y[ind]

M = max(x.max(), y.max())
plt.figure()
plt.plot(x, y, 'o')
plt.plot([0, M], [0, M], 'k')

plt.xlabel("Mean(Author Uncertainty)")
plt.ylabel("Standard Deviation of Measurements")
plt.title("Error Estimates: Density [kg / m^3]")
plt.savefig("./manuscript/figures/error_analysis_density.pdf", bbox_inches=None)

x = data["Relative permittivity at zero frequency_uncertainty_author"]
y = data["Relative permittivity at zero frequency_uncertainty_std"]

ind = x.dropna().index.intersection(y.dropna().index)
x, y = x[ind], y[ind]

M = max(x.max(), y.max())
plt.figure()
plt.plot(x, y, 'o')
plt.plot([0, M], [0, M], 'k')

plt.xlabel("Mean(Author Uncertainty)")
plt.ylabel("Standard Deviation of Measurements")
plt.title("Error Estimates: Relative permittivity at zero frequency")
plt.savefig("./manuscript/figures/error_analysis_dielectric.pdf", bbox_inches=None)


x = data["Mass density, kg/m3_uncertainty_author"]
y = data["Mass density, kg/m3_uncertainty_std"]

plt.figure()
plt.plot(x, 'o', label="mean(author)")
plt.plot(y, 'o', label="std()")

plt.xlabel("Measurment Index")
plt.ylabel("Error Estimate")
plt.title("Error Estimates: Density [kg / m^3]")
plt.legend(loc=0)

plt.savefig("./manuscript/figures/error_analysis_density_index.pdf", bbox_inches=None)



x = data["Relative permittivity at zero frequency_uncertainty_author"]
y = data["Relative permittivity at zero frequency_uncertainty_std"]
plt.figure()
plt.plot(x, 'o', label="mean(author)")
plt.plot(y, 'o', label="std()")

plt.xlabel("Measurement Index")
plt.ylabel("Error Estimate")
plt.title("Error Estimates: Relative permittivity at zero frequency")
plt.legend(loc=0)

plt.savefig("./manuscript/figures/error_analysis_dielectric_index.pdf", bbox_inches=None)
