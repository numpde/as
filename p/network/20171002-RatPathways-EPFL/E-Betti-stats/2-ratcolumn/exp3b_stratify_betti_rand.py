
# RA, 2017-11-17

# Track the Betti numbers in the ratcolumn graph
# as edges are dropped randomly

### IMPORTS -- #

import gc
import pickle
import random
import scipy.io
import numpy     as np
import networkx  as nx
import progressbar, time

from joblib import Parallel, delayed
from topology_localcopy import betti_bin_cpp as betti

### INPUT ---- #

input_file_graph = "../../C-graph1/OUTPUT/UV/column-a-graph.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/exp3b.pkl"

### PARAMS --- #

# Fraction of edges included
FE = np.logspace(-4.2, 0, 31).tolist()

# Desired number of runs for the statistic
max_runs = 20

# Maximal time allowance for each fe (in minutes)
max_time_per_fe = 30

# Number of computing threads to use
num_of_cores = 4

## Use this for testing purposes
#G = nx.gnp_random_graph(40, 0.4, seed=0)

### MEAT ----- #

# See if a graph is already provided
try :
	G
# Otherwise load the ratcolumn graph from disk
except NameError :
	G = pickle.load(open(input_file_graph, "rb"))['G']


# Compute Betti numbers when a fraction of edges are removed randomly
#
# 0 <= fe <= 1 is the remaining fraction of edges
#
def betti_fe(fe) :
	# Number of edges to include
	ne = round(G.number_of_edges() * fe)
	
	# Stratify graph to short edges
	g = nx.Graph()
	g.add_nodes_from(G.nodes())
	g.add_edges_from(random.sample(G.edges(), ne))
	
	# Compute the Betti numbers
	be = betti(nx.find_cliques(g))
	# be[k] is now the k-th Betti number
	
	# Clear
	del g
	gc.collect()
	
	return be

# LL is a list of lists
# Append zeros to each list for uniform length
def padzeros(LL) :
	maxL = max(len(L) for L in LL)
	return [(L + ([0] * (maxL - len(L)))) for L in LL]

# Compute statistics 
def job(fe) :
	BE = []
	
	t0 = time.time()
	for run in range(1, max_runs + 1) :
		BE.append(betti_fe(fe))
		if ((time.time() - t0) > (max_time_per_fe * 60)) : break
	
	BE = np.vstack(padzeros(BE))
	
	(bem, bes) = (np.mean(BE, 0), np.std(BE, 0))
	
	return (bem.tolist(), bes.tolist(), run)

progbar = progressbar.ProgressBar()
BETTI = Parallel(n_jobs=num_of_cores)(delayed(job)(fe) for fe in progbar(FE))

BEM = [bem for (bem, _, _) in BETTI] # Mean
BES = [bes for (_, bes, _) in BETTI] # Std dev
RUN = [run for (_, _, run) in BETTI] # Number of samples

# Collect the results
results = { "FE" : FE, "BEM" : BEM, "BES" : BES, "RUN" : RUN }

# Save to file
pickle.dump(results, open(output_file_stats, "wb"))
