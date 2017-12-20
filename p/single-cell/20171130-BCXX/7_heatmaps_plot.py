
# RA, 2017-12-20

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import numpy as np
import matplotlib.pyplot as plt

## ==================== INPUT :

IFILE = {
	"BC data" : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	"GO -> ENSG" : "OUTPUT/0_e2go/go2e.txt",
}

## =================== OUTPUT :

OFILE = {
}



def heatmap(D) :
	# Source (2017-12-19):
	# https://medium.com/@jerilkuriakose/heat-maps-with-dendrogram-using-python-d112a34e865e

	import numpy as np
	import matplotlib.pyplot as plt
	import scipy.cluster.hierarchy as sch

	# Dendrogram that comes to the left
	fig = plt.figure(figsize=(8,8))
	# Add an axes at position rect [left, bottom, width, height]
	ax1 = fig.add_axes([0.09, 0.1, 0.1, 0.7])
	Y = sch.linkage((D), method='centroid')
	# orientation='left' is reponsible for making the 
	# dendrogram appear to the left
	Z1 = sch.dendrogram(Y, orientation='left')
	ax1.set_xticks([])
	ax1.set_yticks([])

	# top side dendogram
	ax2 = fig.add_axes([0.2, 0.81, 0.7, 0.1])
	Y = sch.linkage(np.transpose(D), method='single')
	Z2 = sch.dendrogram(Y)
	ax2.set_xticks([])
	ax2.set_yticks([])

	# main heat-map
	axmatrix = fig.add_axes([0.2, 0.1, 0.7, 0.7])
	idx1 = Z1['leaves']
	idx2 = Z2['leaves']
	D = D[idx1, :]
	D = D[:, idx2]
	# the actual heat-map
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap="YlGnBu")
	axmatrix.set_xticks([])
	axmatrix.set_yticks([])

	# xticks to the right (x-axis)
	axmatrix.set_xticks(range(D.shape[1]))
	axmatrix.set_xticklabels(idx2, minor=False)
	axmatrix.xaxis.set_label_position('bottom')
	axmatrix.xaxis.tick_bottom()

	plt.xticks(rotation=-90, fontsize=8)

	# xticks to the right (y-axis)
	axmatrix.set_yticks(range(D.shape[0]))
	axmatrix.set_yticklabels(idx1, minor=False)
	axmatrix.yaxis.set_label_position('right')
	axmatrix.yaxis.tick_right()

	# to add the color bar
	axcolor = fig.add_axes([0.01, 0.1, 0.02, 0.7])
	cbar = fig.colorbar(im, cax=axcolor)


# [ LOAD BC DATASET ]

# Load the BC data
BC_data = pickle.load(open(IFILE["BC data"], "rb"))

# Expression matrix
X = BC_data['X']

# ENSG terms
BC_E = BC_data['gene_id']

# Labels for axis/dimension of BC data
(axis_smpl, axis_gene) = (BC_data['axis_smpl'], BC_data['axis_gene'])
	
# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# Batch to indices/header
B2SH = BC_data['B2SH']

# Clusters
S = sorted( tuple(s for (s, h) in SH) for SH in B2SH.values() )


# [ TRANSFORM DATA ]

# Z transform
from scipy import stats
Z = stats.mstats.zscore(X, axis=axis_smpl)


# [ LOAD GO TERMS ]

# Read  GO2E : GO ID --> [ENSG IDs]  from file
GO2E = dict(
	(go_E[0], go_E[1:])
	for go_E in [
		L.rstrip().split('\t') 
		for L in open(IFILE["GO -> ENSG"], 'r')
	]
)


go = 'GO:0006281' # DNA-repair
go = 'GO:0016477' # cell migration
go = 'GO:0001525' # angiogenesis

GO2E = { go : GO2E[go] }

# GO2K is a list of pairs (GO ID, [gene numbers in data])
GO2K = dict( (go, [BC_E.index(e) for e in E if (e in BC_E)]) for (go, E) in GO2E.items() )


# [ SUBSET ]

K = GO2K[go]
X = np.take(X, K, axis=axis_gene)
#print(X.shape)


# [ DIFFERENTIAL EXPRESSION ]

Y = X

# Differential expression within a collections of empirical proba
def KS(P) :
	return np.max([stats.ks_2samp(p, q)[0] for p in P for q in P])

# Reorder X by differentially expressed genes first
DE = []
for k in range(Y.shape[axis_gene]) :
	y = np.take(Y, k, axis=axis_gene)
	DE.append(KS([np.take(y, s, axis=axis_smpl) for s in S]))

I = [n for (n, ks) in sorted(enumerate(DE), key=(lambda x : -x[1]))]

X = X[:, I]


# [ IN-CLUSTER REORDER ]

for s in S :
	Y = X[s, :]
	Y = Y[sorted(range(Y.shape[0]), key=(lambda i : -np.sum(Y[i, :])))]
	X[s, :] = Y



# [ HEATMAP ]

#X = np.take(X, list(range(0, 20)), axis=axis_gene)
#X = np.take(X, list(range(0, 10)), axis=axis_smpl)

Y = np.zeros(X.shape)
#
Y[X != 0] = np.log(X[X != 0])
Y[Y <= 0] = 0
Y = Y / np.max(Y)
#
Y[X == 0] = -1
#
X = Y

#heatmap(X)


R = np.cumsum([0] + [len(s) for s in S])
assert(max(R) == n_samples)

R = R / n_samples * 0.98 + 0.01

fig = plt.figure(figsize=(8,8))

AX = [
	fig.add_axes([a, 0.01, b-a, 0.98])
	for (a, b) in zip(R[:-1], R[1:])
]

for (n, ax) in enumerate(AX) :
	Y = np.take(X, S[n], axis=axis_smpl)
	ax.matshow(np.transpose(Y), aspect='auto', origin='upper', cmap="YlGnBu", vmin=-0.1, vmax=1)
	
	ax.set_xticks([])
	ax.set_yticks([])

plt.show()

