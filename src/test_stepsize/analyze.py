import pandas as pd
import pymbar

table = {}
for timestep in [0.5, 1.0, 1.5, 2.0]:
    max_n = int(5000000 / timestep)
    data = pd.read_csv("./production_%0.2f.log" % timestep)["Density (g/mL)"].values[0:max_n]
    n = len(data)
    g = pymbar.timeseries.statisticalInefficiency(data)
    neff = len(data) / g
    mu = data.mean()
    sigma = data.std()
    stderr = sigma * neff ** -0.5
    table[timestep] = dict(n=n, neff=neff, mu=mu, sigma=sigma, stderr=stderr, timestep=timestep)

table = pd.DataFrame(table).T
table["error"] = table.mu - table.mu[0.5]
table["relerr"] = (table.mu - table.mu[0.5]) / table.mu[0.5]
table
print table.to_latex()
table.to_csv("../../../tables/timestep_dependence.csv")
