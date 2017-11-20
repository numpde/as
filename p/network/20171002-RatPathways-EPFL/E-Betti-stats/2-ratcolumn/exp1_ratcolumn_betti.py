
# RA, 2017-11-06

# Compute the Betti numbers for the ratcolumn graph

### IMPORTS -- #

from topology_localcopy import betti_bin_cpp

import pickle
import gc

### INPUT ---- #
input_file_cliques = "../../C-graph1/OUTPUT/UV/column-b-maxcliques.pkl"

### OUTPUT --- #
pass

### MEAT ----- #

C = pickle.load(open(input_file_cliques, "rb"))['C']

print("Number of maximal cliques:", len(C))

print("Converting C to sorted tuples")
C = [tuple(sorted(mc)) for mc in C]
print("Done")

# Give the garbage collector a hint
gc.collect()

## Test
#import networkx as nx
#G = nx.gnp_random_graph(20, 0.5, seed=0)
#C = nx.find_cliques(G)

print("Betti numbers:", betti_bin_cpp(C, verbose=True))
