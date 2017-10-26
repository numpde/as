
# RA, 20171023

import scipy.io
import pickle
import networkx as nx
import numpy as np

# INPUT
input_filename_graph = "./OUTPUT/UV/column-a-graph.pkl"

# OUTPUT
output_filename_mat  = "./OUTPUT/UV/column-d-laplacian.mat"


# 1. Load the graph

# DEBUG:
#G = nx.complete_graph(10)

try :
	G
except NameError :
	pass
	G = pickle.load(open(input_filename_graph, "rb"))['G']

print("Got G with {} nodes and {} edges".format(G.number_of_nodes(), G.number_of_edges))

# 2. Make the Laplacian

L = nx.laplacian_matrix(G, weight=None)
L = scipy.sparse.csc_matrix(L, dtype=float)
print("Got L of size {}".format(L.shape))

# 3. Save the Laplacian

scipy.io.savemat(output_filename_mat, { 'L' : L })
