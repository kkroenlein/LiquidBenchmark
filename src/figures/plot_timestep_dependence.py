import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("whitegrid")
sns.set(font_scale=1.2)

x = pd.read_csv("./tables/timestep_dependence.csv")


plt.errorbar(x.timestep, x.mu, yerr=x.stderr, fmt='o')
plt.xticks([0.5, 1.0, 2.0])
plt.ylabel("Density [g / cm^3]")
plt.xlabel("Timestep [fs]")
plt.title("Timestep dependence: Density")
plt.savefig("./manuscript/figures/timestep_dependence_density.pdf", bbox_inches="tight")

plt.figure()
reference = x.set_index("timestep").mu[0.5]
plt.errorbar(x.timestep, x.relerr, yerr=x.stderr / reference, fmt='o')
plt.xticks([0.5, 1.0, 2.0])
plt.ylabel("Relative Error [unitless]")
plt.xlabel("Timestep [fs]")
plt.title("Timestep dependence: Relative Error")
plt.savefig("./manuscript/figures/timestep_dependence_relative_error.pdf", bbox_inches="tight")
