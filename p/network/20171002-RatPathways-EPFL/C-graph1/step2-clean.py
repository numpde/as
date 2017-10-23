import scipy.io
import pickle
import networkx as nx
import numpy    as np

def remove_solitary_nodes(G) :
	G.remove_nodes_from([n for (n, d) in G.degree() if (d == 0)])
    
C = pickle.load(open("column-b-cliques.pkl", "rb"))['C']
print("Got C")

(h, _) = np.histogram([len(c) for c in C], bins=range(1, 13))
print("Histogram of clique size:", h)
del h

cc = []
for clique_size in [2, 3] :
	print("Processing clique size:", clique_size)
	
	c0 = []
	for c in C :
		if (len(c) == clique_size) :
			for i in range(0, clique_size) :
				for j in range(i+1, clique_size) :
					c0.append( (c[i], c[j]) )
	
	cc.append(c0)
	del c0

del C


G = pickle.load(open("column-a-graph.pkl", "rb"))['G']
print("Got G")

for c0 in cc : G.remove_edges_from(c0)
print("Removed lower-d clique edges.")

remove_solitary_nodes(G)
print("Removed solitary nodes.")

assert(nx.number_connected_components(G) == 1)

#gg = [g for g in nx.connected_component_subgraphs(G) if (g.number_of_nodes() >= 2)]
#print("Number of nodes in nontrivial connected components:", [g.number_of_nodes() for g in gg])
#print("Minimum edge cut for each:")
#for g in gg : print(nx.minimum_edge_cut(g))

C = list(nx.find_cliques(G))
pickle.dump({'C' : C}, open("column-c-cliques-clean.pkl", "wb"))
print("Done.")
