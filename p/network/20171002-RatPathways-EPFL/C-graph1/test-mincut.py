
# RA, 20171023

import scipy.io
import pickle
import networkx as nx
import numpy    as np

# INPUT
input_filename_graph   = "./OUTPUT/UV/column-a-graph.pkl"



G = pickle.load(open(input_filename_graph, "rb"))['G']
print("Got G")

print("Minimum edge cut:")
print(nx.minimum_edge_cut(G))
