import fcsparser
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import NearestNeighbors as NN

(meta, data) = fcsparser.parse('wishbone_thymus_panel1_rep1.fcs', reformat_meta=True)

print()

# All channel names
channel_names = meta['_channel_names_']

# Channel numbers for channels of interest
CNO = range(0, 10)

X = []
for cno in CNO :
	X.append(data[channel_names[cno]])

X = np.vstack(X).transpose()
#print(X.shape) # X is 250170 x #CNO

nbrs = NN(n_neighbors=2, algorithm='ball_tree').fit(X)
(distances, indices) = nbrs.kneighbors(X)
