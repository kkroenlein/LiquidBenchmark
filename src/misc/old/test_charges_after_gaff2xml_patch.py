import glob, os, chemistry
import networkx
import mdtraj as md

filenames = glob.glob("/home/kyleb/lb_benchmark_openmoltools/tleap/*.prmtop")
#filenames = glob.glob("/home/kyleb/lb_benchmark_openmoltools/tleap/57-55-6_1000_313.2*.prmtop")
filenames = glob.glob("/home/kyleb/lb_benchmark_openmoltools/tleap/121182*.prmtop")
#filenames = glob.glob("/home/kyleb/lb_benchmark_openmoltools/tleap/126492-54-4_1000_315*.prmtop")

for filename in filenames:
    base = os.path.splitext(os.path.split(filename)[-1])[0]
    prmtop_filename = "/home/kyleb/lb_benchmark_openmoltools/tleap/" + base + ".prmtop"
    inpcrd_filename = "/home/kyleb/lb_benchmark_openmoltools/tleap/" + base + ".inpcrd"
    print(prmtop_filename, inpcrd_filename)
    prmtop0 = chemistry.load_file(prmtop_filename)
    t0 = md.load(inpcrd_filename, top=prmtop_filename)
    t0 = t0.atom_slice(t0.top.select("resSeq 0"))
    prmtop_filename = "/home/kyleb/liquid_benchmark_3_14//tleap/" + base + ".prmtop"
    inpcrd_filename = "/home/kyleb/liquid_benchmark_3_14//tleap/" + base + ".inpcrd"
    prmtop1 = chemistry.load_file(prmtop_filename)
    t1 = md.load(inpcrd_filename, top=prmtop_filename)
    t1 = t1.atom_slice(t1.top.select("resSeq 0"))
    top0 = prmtop0.to_dataframe()
    top0 = top0[top0.resid == 0]
    top1 = prmtop1.to_dataframe()
    top1 = top1[top1.resid == 0]
    b0 = t0.top.to_dataframe()[1]
    b1 = t1.top.to_dataframe()[1]
    g0 = networkx.from_edgelist(b0)
    g1 = networkx.from_edgelist(b1)
    for i in range(len(g0.nodes())):
        node = g0.node[i]
        node["charge"] = top0.charge[i]
        node["gafftype"] = top0.type[i]        
    for i in range(len(g1.nodes())):
        node = g1.node[i]
        node["charge"] = top1.charge[i]
        node["gafftype"] = top1.type[i]
    epsilon = 1E-1
    node_match = lambda x, y: (x["gafftype"] == y["gafftype"]) and (abs(x["charge"] - y["charge"]) < epsilon)
    matcher = networkx.isomorphism.GraphMatcher(g0, g1, node_match=node_match)
    matcher.is_isomorphic()
    mapping = matcher.mapping
    mapping
    inv_mapping = {val:key for key, val in mapping.items()}

    top1["newind"] = top1.index.map(lambda x: inv_mapping[x])
    remapped = top1.set_index("newind").sort_index()
    remapped["top0"] = top0["charge"]
    remapped["delta"] = remapped.charge - remapped.top0
    print(remapped.delta.mean(), remapped.delta.std(), remapped.delta.max(), remapped.delta.min())


top0[["name", "type", "charge"]]
remapped[["name", "type", "charge", "top0", "delta"]]
