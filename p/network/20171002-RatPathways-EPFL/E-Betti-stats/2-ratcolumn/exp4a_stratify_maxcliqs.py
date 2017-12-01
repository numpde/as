
# RA, 2017-11-20

# Stats of Betti numbers for subgraphs
# made from a subset of max-cliques

### IMPORTS -- #

import gc
import sys
import pickle
import random
import numpy     as np
import networkx  as nx

from time import time
from datetime import datetime
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress
from topology_localcopy import betti_bin_cpp as betti

### INPUT ---- #

input_file_cliqs = "../../C-graph1/OUTPUT/UV/column-b-maxcliques.pkl"

### OUTPUT --- #

output_file_stats = "./OUTPUT/exp4a_{launched}.pkl"

### PARAMS --- #

# Fraction of cliques included
FC = np.logspace(-3, 0, 31).tolist()

# Desired number of runs for the statistic
n_samples = 100

# Maximal time allowance for each fc (in minutes, after the first sample)
max_time_per_fc = 60

# Number of computing threads to use
max_n_jobs = 8

### MEAT ----- #

# Computation launch time as run identifier
output_file_stats = \
output_file_stats.format(launched = datetime.now().strftime("%Y%m%d-%H%M%S"))

if ("TEST" in sys.argv) :
	# Provide "TEST" on command line for testing purposes
	
	print("# Warning: We are in test mode")
	G = nx.gnp_random_graph(40, 0.8, seed=0)
	C = nx.find_cliques(G)
	n_samples = 4
	max_n_jobs = 2

else :
	# Otherwise load the ratcolumn graph from disk

	#G = pickle.load(open(input_file_graph, "rb"))['G']
	C = pickle.load(open(input_file_cliqs, "rb"))['C']

# Format the max-cliques, commit to RAM
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
	be = betti(c, verbose=("VERBOSE_BETTI" in sys.argv))

	# Clear
	del c
	gc.collect()
	
	return be

def job(fc) :
	return betti_fc(fc)

# Compute statistics over several runs
def stats(fc) :
	
	# Starting time for this round
	t0 = time()

	# Initial number of parallel jobs
	n_jobs = 1
	
	# Result accumulator
	BE = []
	
	while (len(BE) < n_samples) and ((time() - t0) <= (max_time_per_fc * 60)) :
		
		batchsize = min(n_samples - len(BE), n_jobs)
		#print("Doing a batch of {} on {} cores".format(batchsize, n_jobs))
		
		BE.extend(
			Parallel(n_jobs=n_jobs)(
				delayed(job)(fc) for _ in range(batchsize)
			)
		)
		
		# Progressively increase the number of parallel jobs
		n_jobs = min(n_jobs + 1, max_n_jobs)
	
	return BE


BETTI = dict()

for fc in Progress()(FC) :
	BETTI[fc] = stats(fc)

	# Collect the results
	results = { "FC" : FC, "BETTI" : BETTI }

	# Save the (intermediate) results to file
	pickle.dump(results, open(output_file_stats, "wb"))
