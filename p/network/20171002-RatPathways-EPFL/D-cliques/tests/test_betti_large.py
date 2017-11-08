
# RA, 2017-11-06

# See if the boundary operators can be assembled for the ratcolumn graph

import pickle
import gc

C = pickle.load(open("../../C-graph1/OUTPUT/UV/column-b-maxcliques.pkl", "rb"))['C']

print("Number of maximal cliques:", len(C))

import importlib.util as iu
spec = iu.spec_from_file_location("topology", "../topology.py")
topology = iu.module_from_spec(spec)
spec.loader.exec_module(topology)

print("Converting C to sorted tuples")
C = [tuple(sorted(mc)) for mc in C]
print("Done")

# Give the garbage collector a hint
gc.collect()

print("Betti numbers:", topology.betti_bin(C, verbose=True))
