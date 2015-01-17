import pandas as pd

experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]
#experiments = ["Relative permittivity at zero frequency"]

X = pd.read_csv("./tables/full_filtered_data.csv")

data = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].mean().dropna()
counts = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].count().ix[data.index]

uncertainty_std = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].std().ix[data.index]
uncertainty_author = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[[e + "_std" for e in experiments]].mean().ix[data.index]
uncertainty_author.columns = uncertainty_std.columns  # Strip off the _std from column names for now.

# Previously preferred the standard deviation, but now we pick the *larger* error estimate.
#mask = (uncertainty_std.isnull() & (~uncertainty_author.isnull()))
#uncertainty_bestguess = uncertainty_std.copy()
#uncertainty_bestguess[mask] = uncertainty_author[mask]

uncertainty_std.replace(np.nan, -np.inf, inplace=True)
uncertainty_author.replace(np.nan, -np.inf, inplace=True)

uncertainty_bestguess = uncertainty_std.copy()
ind = uncertainty_author > uncertainty_std
uncertainty_bestguess[ind] = uncertainty_author[ind]

uncertainty_std.replace(-np.inf, np.nan, inplace=True)
uncertainty_author.replace(-np.inf, np.nan, inplace=True)
uncertainty_bestguess.replace(-np.inf, np.nan, inplace=True)


for e in experiments:
    data[e + "_uncertainty_std"] = uncertainty_std[e]
    data[e + "_uncertainty_author"] = uncertainty_author[e]
    data[e + "_uncertainty_bestguess"] = uncertainty_bestguess[e]
    data[e + "_counts"] = counts[e]

data.to_csv("./tables/data_with_metadata.csv")
