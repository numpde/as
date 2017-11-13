
# RA, 2017-11-13

# Track the number of cliques in the ratcolumn graph
# as edges are dropped randomly

### IMPORTS -- #

import pickle
import random
import scipy.io
import numpy     as np
import networkx  as nx
import progressbar, time

from collections import Counter
from joblib      import Parallel, delayed

### INPUT ---- #

input_file_graph = "../C-graph1/OUTPUT/UV/column-a-graph.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/column-stratify-stats-2b-rand.pkl"

### PARAMS --- #

# Fraction of edges included
FE = [i/100 for i in range(0, 101)]

# Maximal time allowance for each fe (in seconds)
max_time_per_fe = 1000

# Maximal number of runs for the statistic
max_runs = 100

# Number of computing cores to use
# Each thread requires about 30-50 GB of RAM
num_of_cores = 4

## Use this for testing purposes
#G = nx.gnp_random_graph(100, 0.1, seed=0)

### MEAT ----- #

# See if a graph is already provided
try :
	G
# Otherwise load the ratcolumn graph from disk
except NameError :
	G = pickle.load(open(input_file_graph, "rb"))['G']

# Count cliques when a fraction of edges are removed randomly
#
# 0 <= fe <= 1 is the remaining fraction of edges
#
def cliques(fe) :
	# Number of edges to include
	ne = round(G.number_of_edges() * fe)
	# Stratify graph to short edges
	g = nx.Graph()
	g.add_nodes_from(G.nodes())
	g.add_edges_from(random.sample(G.edges(), ne))
	# Count k-cliques
	C = nx.enumerate_all_cliques(g)
	nc = dict(Counter(len(c) for c in C))
	# Convert nc to a list 
	nc = [nc.get(k, 0) for k in range(0, 1 + max(nc.keys()))]
	# nc[k] is now the number of k-cliques
	return nc

# LL is a list of lists
# Append zeros to each list for uniform length
def padzeros(LL) :
	maxL = max(len(L) for L in LL)
	return [(L + ([0] * (maxL - len(L)))) for L in LL]

# Compute statistics 
def job(fe) :
	NC = []
	
	t0 = time.time()
	for run in range(0, max_runs) :
		NC.append(cliques(fe))
		if ((time.time() - t0) > max_time_per_fe) : break
	
	NC = np.vstack(padzeros(NC))
	
	(ncm, ncs) = (np.mean(NC, 0), np.std(NC, 0))
	
	return (ncm.tolist(), ncs.tolist(), run)

progbar = progressbar.ProgressBar()
CLIQUES = Parallel(n_jobs=num_of_cores)(delayed(job)(fe) for fe in progbar(FE))

NCM = [ncm for (ncm, _, _) in CLIQUES] # Mean
NCS = [ncs for (_, ncs, _) in CLIQUES] # Std dev
RUN = [run for (_, _, run) in CLIQUES] # Number of samples

# Collect the results
results = { "FE" : FE, "NCM" : NCM, "NCS" : NCS, "RUN" : RUN }

# Save to file
pickle.dump(results, open(output_file_stats, "wb"))
