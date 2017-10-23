import scipy.io
import pickle
import networkx as nx
import numpy    as np


print("1. Loading adjacency matrix.")
M = scipy.io.loadmat("../A-h5-to-txt/pathways_mc0_Column.h5.mat")["M"]
print("Done.")



print("--------------------------------")



print("2. Constructing the connectivity graph.")

G = nx.Graph()
G.add_edges_from(zip(*np.nonzero(M)))
pickle.dump({'G' : G}, open("column-a-graph.pkl", "wb"))
print("Done.")



print("--------------------------------")



N = 100
print("3. Extracting a subgraph on {} random nodes.".format(N))

G1 = G.subgraph(np.random.choice(G.nodes(), N))
print("Done.")


#print("HACK! Replacing graph by subgraph"); G = G1



print("--------------------------------")



print("4. Looking for maximal cliques in the subgraph.")

C = list(nx.find_cliques(G1))
print("Done.")

cc = [len(c) for c in C]
(h, _) = np.histogram(cc, bins=range(1, 10))
print("Found: {} cliques.".format(len(cc)))
print("Histogram of clique size:", h)



print("--------------------------------")



print("5. Looking for maximal cliques in the whole graph.")

C = list(nx.find_cliques(G))
pickle.dump({'C' : C}, open("column-b-cliques.pkl", "wb"))
print("Done.")

cc = [len(c) for c in C]
(h, _) = np.histogram(cc, bins=range(1, 13))
print("Found: {} cliques.".format(len(cc)))
print("Histogram of clique size:", h)



