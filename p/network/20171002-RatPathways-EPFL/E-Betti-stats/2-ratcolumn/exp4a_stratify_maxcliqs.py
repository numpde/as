
# RA, 2017-11-20

# Stats of Betti numbers for subgraphs
# made from a subset of max-cliques

### IMPORTS -- #

import gc
import pickle
import random
import numpy     as np
import networkx  as nx

from time import time
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress
from topology_localcopy import betti_bin_cpp as betti

### INPUT ---- #

input_file_cliqs = "../../C-graph1/OUTPUT/UV/column-b-maxcliques.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/exp4a.pkl"

### PARAMS --- #

# Fraction of cliques included
FC = np.logspace(-3, 0, 31).tolist()

# Desired number of runs for the statistic
max_runs = 100

# Maximal time allowance for each fc (in minutes, after first sample)
max_time_per_fc = 60

# Number of computing threads to use
num_of_cores = 8

## Use this for testing purposes
#G = nx.gnp_random_graph(40, 0.8, seed=0)
#C = nx.find_cliques(G)

### MEAT ----- #

# See if a graph is already provided
try :
	G
# Otherwise load the ratcolumn graph from disk
except NameError :
	C = pickle.load(open(input_file_cliqs, "rb"))['C']

# Get max-cliques
C = [tuple(sorted(mc)) for mc in C]

# Compute Betti numbers for a graph on a subset of cliques
#
# 0 <= fc <= 1 is the remaining fraction of cliques
#
def betti_fc(fc) :
	# Number of max-cliques to include
	nc = round(len(C) * fc)
	
	# Choose nc max-cliques
	c = random.sample(C, k=nc)
	
	# Compute the Betti numbers
	# be[k] is the k-th Betti number
	be = betti(c, verbose=(num_of_cores==1))

	# Clear
	del c
	gc.collect()
	
	return be

# LL is a list of lists
# Append zeros to each list for uniform length
def padzeros(LL) :
	maxL = max(len(L) for L in LL)
	return [(L + ([0] * (maxL - len(L)))) for L in LL]

def job(fc, t0) :
	if ((time() - t0) > (max_time_per_fc * 60)) : return None
	
	return betti_fc(fc)

# Compute statistics over several runs
def stats(fc) :
	BE = Parallel(n_jobs=num_of_cores)(delayed(job)(fc, time()) for _ in range(max_runs))
	
	BE = [be for be in BE if (be is not None)]
	
	run = len(BE) # Number of samples
	assert(run), "Failed to get any samples."
	
	BE = np.vstack(padzeros(BE))
	
	(bem, bes) = (np.mean(BE, 0), np.std(BE, 0))
	
	return (bem.tolist(), bes.tolist(), run)

BETTI = [stats(fc) for fc in Progress()(FC)]

BEM = [bem for (bem, _, _) in BETTI] # Mean
BES = [bes for (_, bes, _) in BETTI] # Std dev
RUN = [run for (_, _, run) in BETTI] # Number of samples

# Collect the results
results = { "FC" : FC, "BEM" : BEM, "BES" : BES, "RUN" : RUN }

# Save to file
pickle.dump(results, open(output_file_stats, "wb"))
