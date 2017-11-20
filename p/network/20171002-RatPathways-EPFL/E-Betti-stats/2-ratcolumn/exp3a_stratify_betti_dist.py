
# RA, 2017-11-16

# Track the Betti numbers in the ratcolumn graph
# as edges between distant nodes are dropped

### IMPORTS -- #

import pickle
import scipy.io
import progressbar
import numpy    as np
import networkx as nx

from joblib import Parallel, delayed
from topology_localcopy import betti_bin_cpp as betti

### INPUT ---- #

input_file_xyz = "../../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat"
input_file_graph = "../../C-graph1/OUTPUT/UV/column-a-graph.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/exp3a.pkl"

### PARAMS --- #

# Fraction of maximal edge length
P = np.logspace(-3, 0, 31).tolist()

# Number of computing threads to use in parallel
num_of_cores = 10

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

# Compute Betti numbers when long edges are removed
#
# The parameter 0 <= p <= 1 is the cut-off
# of edge length wrt the longest edge
#
def betti_p(p) :
	# Maximal allowed edge length
	d = maxd * p
	# Stratify graph to short edges
	g = nx.Graph()
	g.add_nodes_from(G.nodes())
	g.add_edges_from((a, b) for (a, b, data) in G.edges(data=True) if (data['d'] <= d))
	# Fraction of edges included
	fe = g.number_of_edges() / G.number_of_edges()
	# Compute the Betti numbers
	be = betti(nx.find_cliques(g))
	# be[k] is now the k-th Betti number
	return (fe, be)

def job(p) :
	return betti_p(p)

progbar = progressbar.ProgressBar()
BETTI_P = Parallel(n_jobs=num_of_cores)(delayed(job)(p) for p in progbar(P))

FE = [fe for (fe, _) in BETTI_P]
BE = [be for (_, be) in BETTI_P]

# Collect the results
results = { "maxd" : maxd, "P" : P, "FE" : FE, "BE" : BE }

# Save to file
pickle.dump(results, open(output_file_stats, "wb"))
