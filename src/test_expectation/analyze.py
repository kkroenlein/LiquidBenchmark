import pandas as pd
import pymbar

data1 = pd.read_csv("./production_langevin1fs.log")["Density (g/mL)"].values
data2 = pd.read_csv("./production_langevin.log")["Density (g/mL)"].values

g1 = pymbar.timeseries.statisticalInefficiency(data1)
neff1 = len(data1) / g1
g2 = pymbar.timeseries.statisticalInefficiency(data2)
neff2 = len(data2) / g2

mu1 = data1.mean()
mu2 = data2.mean()

sigma1 = data1.std()
sigma2 = data2.std()

err1 = sigma1 * neff1 ** -0.5
err2 = sigma2 * neff2 ** -0.5

(mu1 - mu2)
(mu1 - mu2) / err1
