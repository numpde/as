
# RA, 2017-11-11

# Track the number of cliques in the ratcolumn graph
# as edges between distant nodes are dropped

### IMPORTS -- #

import pickle
import scipy.io
import progressbar
import numpy    as np
import networkx as nx

from joblib      import Parallel, delayed
from collections import Counter

### INPUT ---- #

input_file_xyz = "../../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat"
input_file_graph = "../../C-graph1/OUTPUT/UV/column-a-graph.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/column-stratify-stats-2a-dist.pkl"

### PARAMS --- #

# Fraction of maximal edge length
P = np.logspace(-3, 0, 101).tolist()

# Number of computing cores to use
# Each thread requires 30-50 GB of RAM
num_of_cores = 4

## Use this for testing purposes
#G = nx.gnp_random_graph(100, 0.4, seed=0)
#for (a, b, ed) in G.edges(data=True) : ed['d'] = np.random.rand()

### MEAT ----- #

# See if a graph is already provided
try :
	G
# Otherwise load the ratcolumn graph from disk
except NameError :
	G = pickle.load(open(input_file_graph, "rb"))['G']

	# Nodes locations in 3D
	X = scipy.io.loadmat(input_file_xyz)["XYZ"]

	# Compute edge lengths; use euclidean distance
	for (a, b, ed) in G.edges(data=True) :
		ed['d'] = np.linalg.norm(X[a, :] - X[b, :], 2)

# Maximal edge length
maxd = max(ed['d'] for (a, b, ed) in G.edges(data=True))

# Count cliques when long edges are removed
#
# The parameter 0 <= p <= 1 is the cut-off
# of edge length wrt the longest edge
#
def cliques(p) :
	# Maximal allowed edge length
	d = maxd * p
	# Stratify graph to short edges
	g = nx.Graph()
	g.add_nodes_from(G.nodes())
	g.add_edges_from((a, b) for (a, b, data) in G.edges(data=True) if (data['d'] <= d))
	# Fraction of edges included
	fe = g.number_of_edges() / G.number_of_edges()
	# Count k-cliques
	C = nx.enumerate_all_cliques(g)
	nc = dict(Counter(len(c) for c in C))
	# Convert nc to a list 
	nc = [nc.get(k, 0) for k in range(0, 1 + max(nc.keys()))]
	# nc[k] is now the number of k-cliques
	return (fe, nc)

def job(p) :
	return cliques(p)

progbar = progressbar.ProgressBar()
CLIQUES = Parallel(n_jobs=num_of_cores)(delayed(job)(p) for p in progbar(P))

FE = [fe for (fe, _) in CLIQUES]
NC = [nc for (_, nc) in CLIQUES]

# Collect the results
results = { "maxd" : maxd, "P" : P, "FE" : FE, "NC" : NC }

# Save to file
pickle.dump(results, open(output_file_stats, "wb"))
