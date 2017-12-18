
# RA, 2017-11-20

# Properties of the max-cliques hypergraph

### IMPORTS -- #

import numpy     as np
import networkx  as nx

### INPUT ---- #

input_file_graph = "../../C-graph1/OUTPUT/UV/column-a-graph.pkl"

### OUTPUT --- #

#output_file_stats = "./OUTPUT/exp1.pkl"

### PARAMS --- #

# Use this for testing purposes
G = nx.gnp_random_graph(40, 0.4, seed=0)

### MEAT ----- #

# See if a graph is already provided
try :
	G
# Otherwise load the ratcolumn graph from disk
except NameError :
	G = pickle.load(open(input_file_graph, "rb"))['G']


## https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.clique.make_clique_bipartite.html
#B = nx.make_clique_bipartite(G)

# https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.clique.make_max_clique_graph.html
H = nx.make_max_clique_graph(G)

# https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.clique.find_cliques.html
C = list(nx.find_cliques(G))

# We assume that H and C are compatible:
#     a and b is are connected nodes in H 
#     <==>
#     C[a] and C[b] have common vertices.
#
# Let the weight 'w' of a hyperedge be
# the number of common vertices
#
# Let the weight 'w' of a hypernode be
# the number of vertices in the max-clique

for (a, b, d) in H.edges(data=True) :
	d['w'] = len(set(C[a]).intersection(set(C[b])))

for (a, d) in H.nodes(data=True) :
	d['w'] = len(C[a])

# Largest max-clique size
K = max(len(C[a]) for a in H.nodes())

# Histogram of max-clique sizes
count = np.zeros(K+1, dtype='int')
# Fill the Histogram
for mc in C : count[len(mc)] += 1

for (k, c) in list(enumerate(count))[1:] :
	print("Number of maximal {}-cliques:".format(k), c)

# 
degrees = np.zeros( (K+1, K+1), dtype='int' )

for (a, b) in H.edges() :
	(m, n) = (len(C[a]), len(C[b]))
	
	degrees[m, n] += 1
	degrees[n, m] += 1

print(degrees)
