
# RA, 20171023

import scipy.io
import pickle
import networkx as nx

# INPUT
input_filename_graph   = "./OUTPUT/UV/column-a-graph.pkl"



G = pickle.load(open(input_filename_graph, "rb"))['G']
print("Got G")


L = nx.laplacian_matrix(G, weight=None)
print("Got L")


print("Laplacian synopsis:", L.shape)

