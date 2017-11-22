
# RA, 2017-10-23

### IMPORTS -- #

import gc
import scipy.io
import pickle
import networkx as nx
import numpy    as np

### INPUT ---- #

input_filename = "../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat"

### OUTPUT --- #

output_filename_graph      = "./OUTPUT/UV/column-a-graph.pkl"
output_filename_maxcliques = "./OUTPUT/UV/column-b-maxcliques.pkl"
#output_filename_allcliques = "./OUTPUT/UV/column-b-allcliques.pkl"

### MEAT ----- #

print("1. Loading adjacency matrix.")
M = scipy.io.loadmat(input_filename)["M"]
print("Done.")


print("--------------------------------")


print("2. Constructing the connectivity graph.")

G = nx.Graph()
G.add_edges_from(zip(*np.nonzero(M)))
pickle.dump({'G' : G}, open(output_filename_graph, "wb"))
print("Done.")


print("--------------------------------")


N = 100
print("3. Extracting a subgraph on {} random nodes (just a test).".format(N))

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


try :
	# If variable exists ...
	output_filename_maxcliques
	
	# ... continue:

	print("5. Looking for maximal cliques in the whole graph.")

	C = list(nx.find_cliques(G))
	pickle.dump({'C' : C}, open(output_filename_maxcliques, "wb"))
	print("Done.")

	cc = [len(c) for c in C]
	(h, _) = np.histogram(cc, bins=range(1, 13))
	print("Found: {} max cliques.".format(len(cc)))
	print("Histogram of max clique size:", h)


except NameError :
	pass


print("--------------------------------")


try :
	# If variable exists ...
	output_filename_allcliques
	
	# ... continue:

	print("6. Looking for all cliques in the whole graph.")
	
	del C
	gc.collect()
	
	S = []
	for c in nx.enumerate_all_cliques(G) :
		while (len(S) < len(c)) : S.append([])
		S[len(c) - 1].append(tuple(sorted(c)))
	
	pickle.dump({'S' : S}, open(output_filename_allcliques, "wb"))
	print("Done.")
	
	cc = [len(c) for c in C]
	(h, _) = np.histogram(cc, bins=range(1, 13))
	print("Found: {} cliques.".format(len(cc)))
	print("Histogram of clique size:", h)

except NameError :
	pass
