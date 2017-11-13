
# RA, 2017-11-04

# Statistics for Betti numbers of the random G_{n,p} graph

### IMPORTS -- #

import networkx as nx
import numpy    as np

import progressbar, time

from topology_localcopy import betti_bin_cpp as BETTI

# https://blog.dominodatalab.com/simple-parallelization/
from joblib import Parallel, delayed

### INPUT ---- #
pass

### PARAMS --- #

# The highest Betti number rank that might occur
max_expected_betti = 10
# Number of nodes in the random graph
nn = [10, 20, 30, 40]
# Maximal number of runs for the statistic
max_runs = 1000
# Maximal time allowance for each p (in seconds)
max_time_per_p = 1000
# Range of edge proba in G_{n,p}
P = [p/100 for p in range(1, 77)]
# Number of computing cores to use
num_of_cores = 4

### OUTPUT --- #

output_file_stats = "./OUTPUT/gnpbetti-out_n={}.pkl"

### MEAT ----- #

# Computes average and std dev of Betti numbers
# for the G_{n,p} random graph over (runs) runs
def betti_av(n, p, runs) :
	#print("n: {}, p: {}".format(n, p))
	
	B = np.zeros( (0, 1 + max_expected_betti) )
	
	t0 = time.time()
	for i in range(0, runs) :
		G = nx.fast_gnp_random_graph(n, p, seed=i)
		b = BETTI(nx.find_cliques(G))
		# Append to B, padding with zeros
		B = np.vstack([B, b + ([0] * (B.shape[1] - len(b)))])
		
		# Terminate if it takes too long
		if ((time.time() - t0) > max_time_per_p) : break

	return (np.mean(B, 0), np.std(B, 0), i)

# Pickable job wrapper
def job(p) :
	return betti_av(n, p, max_runs)

def stats(n) :
	print("Collecting stats for n = {}".format(n))
	
	# Compute the Betti numbers as (average, stddev) for each p
	bar = progressbar.ProgressBar()
	# https://blog.dominodatalab.com/simple-parallelization/
	MS = Parallel(n_jobs=num_of_cores)(delayed(job)(p) for p in bar(P))

	# Separate the average
	M = [B[0] for B in MS]
	# Separate the standard deviation
	S = [B[1] for B in MS]
	# Separate the number of runs taken
	R = [B[2] for B in MS]

	# Write output
	import pickle
	data = {'M' : M, 'S' : S, 'n' : n, 'P' : P, 'runs' : R}
	filename = output_file_stats.format(n)
	print("Writing results to:", filename)
	pickle.dump(data, open(filename, "wb"))

if __name__ == "__main__":
	
	for n in nn : stats(n)

