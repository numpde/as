
# RA, 2017-11-04

import topology
import networkx as nx

G0 = nx.fast_gnp_random_graph(20, 0.5, seed=255)

G1 = nx.make_max_clique_graph(G0)


print(sorted([d for (n, d) in G0.degree()], reverse=True))
print(sorted([d for (n, d) in G1.degree()], reverse=True))
