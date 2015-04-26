"""
Test building mol2 files using latest gaff2xml build.
"""
import pandas as pd
import chemistry
from simtk import unit as u
from liquid_tools import AmberMixtureSystem
from density_simulation_parameters import MOLECULES_PER_BOX, DATA_PATH
import gaff2xml

data = pd.read_csv("../../tables/data_dielectric.csv")

for k0, k1, components, smiles, cas, temperature, pressure, density, perm in data.itertuples():
    print(k0, k1, components, smiles, cas, temperature, pressure, density)
    mol2_filename = "./tmpmol2/%s.mol2" % cas
    frcmod_filename = "./tmpmol2/%s.frcmod" % cas
    gaff2xml.openeye.smiles_to_antechamber(smiles, mol2_filename, frcmod_filename)



PATH0 = "/home/kyleb/liquid_benchmark_4_8/"
PATH1 = "/home/kyleb/liquid_benchmark_3_14/"
PATH2 = "./tmpmol2/"

for k0, k1, components, smiles, cas, temperature, pressure, density, perm in data.itertuples():
    print(k0, k1, components, smiles, cas, temperature, pressure, density)

    mol2_filename = PATH0 + "/monomers/%s.mol2" % cas
    prmtop0 = chemistry.load_file(mol2_filename)
    charges0 = np.array([a.charge for a in prmtop0.atoms])

    mol2_filename = PATH1 + "/monomers/%s.mol2" % cas
    prmtop1 = chemistry.load_file(mol2_filename)
    charges1 = np.array([a.charge for a in prmtop1.atoms])
    
    mol2_filename = PATH2 + "/%s.mol2" % cas
    prmtop2 = chemistry.load_file(mol2_filename)
    charges2 = np.array([a.charge for a in prmtop2.atoms])
    print(cas)
    print(np.std((charges0.min(), charges1.min(), charges2.min())))
    print(charges0.min(), charges1.min(), charges2.min())
