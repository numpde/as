
# RA, 2017-11-15

import importlib.util as iu
spec = iu.spec_from_file_location("topology", "../topology.py")
topology = iu.module_from_spec(spec)
spec.loader.exec_module(topology)

def betti(C) :
	return topology.betti_bin_cpp(C, worker="../cpp/UV/rank")
	#return topology.betti_bin(C)

import networkx as nx
import numpy as np
import itertools

G = nx.fast_gnp_random_graph(15, 0.6, seed=0)

C = list(nx.find_cliques(G))
H = nx.make_max_clique_graph(G)
assert(H.number_of_nodes() == len(C))

for (a, b) in H.edges() :
	assert(len(set(C[a]).intersection(C[b])))

print("Full:")
print("G:", betti(nx.find_cliques(G)))
print("H:", betti(nx.find_cliques(H)))

H.remove_nodes_from([n for n in H.nodes() if (len(C[n]) <= 3)])
print("After removal:")
print("G:", betti(nx.find_cliques(G)))
print("H:", betti(nx.find_cliques(H)))
