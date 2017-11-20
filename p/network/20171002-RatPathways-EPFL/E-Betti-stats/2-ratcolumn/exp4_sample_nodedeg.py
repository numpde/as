
# RA, 2017-11-20

# Stats of Betti numbers for random graphs
# whose node-degrees approximate those
# of the rat column graph

### IMPORTS -- #

import gc
import pickle
import random
import numpy     as np
import networkx  as nx
import progressbar, time

from joblib import Parallel, delayed
from topology_localcopy import betti_bin_cpp as betti

### INPUT ---- #

input_file_graph = "../../C-graph1/OUTPUT/UV/column-a-graph.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/exp4.pkl"

### PARAMS --- #

# Fraction of nodes included
FN = np.logspace(-1, 0, 11).tolist()

# Desired number of runs for the statistic
max_runs = 100

# Maximal time allowance for each fn (in minutes, after first sample)
max_time_per_fn = 60

# Number of computing threads to use
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

# Node degrees
degrees = list(d for (_, d) in G.degree())

# Compute Betti numbers for a graph on (fn * n) nodes
# whose node-degrees approximate those of the parent graph G
#
# 0 <= fn <= 1 is the remaining fraction of nodes
#
def betti_fn(fn) :
	# Number of edges to include
	nn = round(G.number_of_nodes() * fn)
	
	for _ in range(0, 1000) :
		# Node-degrees sample
		deg = random.sample(degrees, k=nn)
		
		# Can a graph with these node-degrees be constructed?
		if nx.is_graphical(deg) : break
	
		# If not, the sequence is invalid
		deg = None
	
	# Could not find a valid node-degree sample
	if (deg is None) : return None
	
	# Smaller graph, with comparable node degree distribution
	g = nx.havel_hakimi_graph(deg)
	
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
def job(fn) :
	BE = []
	
	t0 = time.time()
	for run in range(1, max_runs + 1) :
		BE.append(betti_fn(fn))
		if ((time.time() - t0) > (max_time_per_fn * 60)) : break
	
	BE = [be for be in BE if (be is not None)]
	
	run = len(BE)
	BE = np.vstack(padzeros(BE))
	
	(bem, bes) = (np.mean(BE, 0), np.std(BE, 0))
	
	return (bem.tolist(), bes.tolist(), run)

progbar = progressbar.ProgressBar()
BETTI = Parallel(n_jobs=num_of_cores)(delayed(job)(fn) for fn in progbar(FN))

BEM = [bem for (bem, _, _) in BETTI] # Mean
BES = [bes for (_, bes, _) in BETTI] # Std dev
RUN = [run for (_, _, run) in BETTI] # Number of samples

# Collect the results
results = { "FN" : FN, "BEM" : BEM, "BES" : BES, "RUN" : RUN }

# Save to file
pickle.dump(results, open(output_file_stats, "wb"))
