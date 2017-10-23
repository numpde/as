import scipy.io
import pickle
import networkx as nx
import numpy    as np


G = pickle.load(open("column-a-graph.pkl", "rb"))['G']
print("Got G")

print("Minimum edge cut:")
print(nx.minimum_edge_cut(G))
