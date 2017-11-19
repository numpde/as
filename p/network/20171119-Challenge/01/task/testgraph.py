
def make() : 
	import networkx as nx
	# Erdos-Renyi random graph
	n = 40
	p = 0.7
	G = nx.gnp_random_graph(n, p, seed=0)
	return G

