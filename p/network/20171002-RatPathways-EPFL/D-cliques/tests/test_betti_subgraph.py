
# RA, 2017-11-04

import topology

import random
import networkx as nx

G = nx.fast_gnp_random_graph(20, 0.5, seed=255)
C = list(nx.find_cliques(G))
B = topology.betti(C)
print("Betti numbers:", B)

H = nx.make_max_clique_graph(G)
#print(H.nodes(), H.edges())

#for (a, b) in H.edges() : print(list(set(C[a]) & set(C[b])))

C3 = [c for c in C if (len(c) == 4)]
print(len(C3))
print("Betti numbers of C3:", topology.betti(C3))
for _ in range(0, 1000) :
	n = 10
	#subC = [C[i] for i in random.sample(H.nodes(), n)]
	#subC = random.sample(C3, n)
	subC = nx.find_cliques(nx.subgraph(G, random.sample(G.nodes(), n)))
	#print(subC)
	B = topology.betti(subC)
	print("Betti numbers of subgraph:", B)

