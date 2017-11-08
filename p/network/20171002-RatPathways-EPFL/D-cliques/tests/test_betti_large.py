
# RA, 2017-11-06

# See if the boundary operators can be assembled for the ratcolumn graph

import pickle

C = pickle.load(open("../../../C-graph1/OUTPUT/UV/column-b-cliques.pkl", "rb"))['C']

print("Number of maximal cliques:", len(C))


import importlib.util as iu
spec = iu.spec_from_file_location("topology", "../../topology.py")
topology = iu.module_from_spec(spec)
spec.loader.exec_module(topology)

print("Betti numbers:", topology.betti(C, verbose=True))
