
# RA, 2017-11-04

# Statistics for Betti numbers of the random G_{n,p} graph

import importlib.util as iu
spec = iu.spec_from_file_location("topology", "../topology.py")
topology = iu.module_from_spec(spec)
spec.loader.exec_module(topology)

import networkx as nx
import numpy as np

def betti_av(n, p, runs) :
	print("n: {}, p: {}".format(n, p))
	B = np.zeros( (0, 10) )
	
	for _ in range(0, runs) :
		G = nx.fast_gnp_random_graph(n, p)
		b = topology.betti(nx.find_cliques(G))
		# Append to B, padding with zeros
		B = np.vstack([B, b + ([0] * (B.shape[1] - len(b)))])

	return (np.mean(B, 0), np.std(B, 0))


# Number of nodes
n = 20
# Number of runs for the statistic
runs = 10000

P = [p/100 for p in range(1, 99)]

MS = [betti_av(n, p, runs) for p in P]

# Average
M = [B[0] for B in MS]
# Standard deviation
S = [B[1] for B in MS]


import pickle
data = {'M' : M, 'S' : S, 'n' : n, 'P' : P, 'runs' : runs}
pickle.dump(data, open("test_gnpbetti.py-out.pkl", "wb"))

#import matplotlib.pyplot as plt
#plt.plot(P, M)
#plt.legend(["b{}".format(n) for n in range(0, 10)])
#plt.show()
