
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

from scipy import stats
from joblib import Parallel, delayed
from itertools import chain
from progressbar import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	'BC data' : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'gene-ks' : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
}

## =================== OUTPUT :

OFILE = {
	'heatmap' : "OUTPUT/7_heatmaps/heatmap_{go}.{ext}"
}


## ====================== (!) :


# [ LOAD BC DATASET ]

# Load the BC data
BC_data = pickle.load(open(IFILE['BC data'], 'rb'))

# Expression matrix
X = BC_data['X']

# Labels for axis/dimension of BC data
(axis_smpl, axis_gene) = (BC_data['axis_smpl'], BC_data['axis_gene'])
	
# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# ENSG IDs
BC_E = BC_data['gene_id']
assert(len(BC_E) == X.shape[axis_gene]), "Inconsistent gene info"

# Batch to indices/header
B2SH = BC_data['B2SH']

# Clusters
G = sorted(list(B2SH.keys()))
S = [ np.array([s for (s, h) in B2SH[g]]) for g in G ]


# [ LOAD KS DATA ]

KS_data = pickle.load(open(IFILE['gene-ks'], 'rb'))
#print(KS.keys())

E2KS = KS_data['E2KS']
KS_meta = KS_data['KS_meta']

## [ TRANSFORM DATA ]

## Z transform
#from scipy import stats
#Z = stats.mstats.zscore(X, axis=axis_smpl)


# [ LOAD GO TERMS ]

# Read  GO2E : GO ID --> [ENSG IDs]  from file
GO2E = {
	go_E[0] : go_E[1:]
	for go_E in [
		L.rstrip().split('\t') 
		for L in open(IFILE['GO=>ENSG'], 'r')
	]
}


## [ PLOT EXPRESSION VS DE ]

#E2X = dict(zip(BC_E, np.moveaxis(X, axis_gene, 0)))

##print(len(E2KS), len(E2X), len(set(E2KS.keys()) & set(E2X.keys())), len(BC_E))

#ks_op = 1

#ks = [ E2KS[e][ks_op]  for e in BC_E ]
#ex = [ np.mean(E2X[e]) for e in BC_E ]

#plt.clf()
#plt.loglog(ex, ks, 'b.', markersize=0.5)
#plt.xlabel('Average expression')
#plt.ylabel('Differential expression (KS {})'.format(KS_meta[ks_op]))
#plt.show()


# [ PLOT EXPRESSION VS DE / GO ]

E2X = dict(zip(BC_E, np.moveaxis(X, axis_gene, 0)))

ks_op = 1

plt.clf()

E = BC_E
ks = [ E2KS[e][ks_op]  for e in E ]
ex = [ np.std(E2X[e]) for e in E ]
plt.loglog(ex, ks, 'g.', markersize=1)

go = 'GO:0006281' # DNA-repair
E = set(GO2E[go]) & set(BC_E)
ks = [ E2KS[e][ks_op]  for e in E ]
ex = [ np.std(E2X[e]) for e in E ]
plt.loglog(ex, ks, 'rd', markersize=3)

go = 'GO:0001525' # angiogenesis
E = set(GO2E[go]) & set(BC_E)
ks = [ E2KS[e][ks_op]  for e in E ]
ex = [ np.std(E2X[e]) for e in E ]
plt.loglog(ex, ks, 'bo', markersize=3)

go = 'GO:0004984' # olfactory receptor activity
E = set(GO2E[go]) & set(BC_E)
ks = [ E2KS[e][ks_op]  for e in E ]
ex = [ np.std(E2X[e]) for e in E ]
plt.loglog(ex, ks, 'kx', markersize=3)

plt.xlabel('Expression std dev')
plt.ylabel('Differential expression (KS {})'.format(KS_meta[ks_op]))
plt.show()


# [ SELECT GO CATEGORY ]

go = 'GO:0001525' # angiogenesis
go = 'GO:0006281' # DNA-repair
go = 'GO:0016477' # cell migration

GO2E = { go : GO2E[go] }

# GO2K : GO ID --> [gene numbers in data]
GO2K = {
	go : [BC_E.index(e) for e in (set(E) & set(BC_E))]
	for (go, E) in GO2E.items() 
}


# [ SUBSET ]

K = GO2K[go]
X = np.take(X, K, axis=axis_gene)
#print(X.shape)


# [ DIFFERENTIAL EXPRESSION ]

Y = X

# Reorder X by differentially expressed genes first
DE = [
	KS([y[s] for s in S])
	for y in np.moveaxis(Y, axis_gene, 0)
]

I = [n for (n, ks) in sorted(enumerate(DE), key=(lambda x : -x[1]))]

X = X[:, I]


# [ IN-CLUSTER REORDER ]

for s in S :
	Y = X[s, :]
	Y = Y[sorted(range(Y.shape[0]), key=(lambda i : -np.sum(Y[i, :])))]
	X[s, :] = Y

# The in-cluster reshuffling invalidates the sample labels
if 'header' in locals() : del header


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

#heatmap(X)


R = [len(s) for s in S]
R = list(chain.from_iterable([[2, b] for b in R]))
R = np.cumsum(R)
R = R / max(R)
R = R * 0.98 + 0.01

phi = (1 + math.sqrt(5)) / 2
fig = plt.figure(figsize=(8*phi, 8))

AX = [
	fig.add_axes([a, 0.01, b-a, 0.98])
	for (a, b) in zip(R[0::2], R[1::2])
]

for (n, ax) in enumerate(AX) :
	y = np.take(Y, S[n], axis=axis_smpl)
	ax.matshow(np.transpose(y), aspect='auto', origin='upper', cmap="YlGnBu", vmin=-0.1, vmax=1)
	
	ax.set_xticks([])
	ax.set_yticks([])
	
	ax.text(0.95, 0.005, G[n], rotation=90, transform=ax.transAxes, horizontalalignment='right', verticalalignment='bottom')

plt.savefig(OFILE['heatmap'].format(go=go.replace(':', '-'), ext='eps'))
plt.show()


plt.clf()
	


## [ PLOT EXPRESSION PDF ]

#plt.clf()
#for i in range(X.shape[axis_gene]) :
	#plt.clf()
	#x = np.take(X, i, axis=axis_gene)
	#for (n, s) in enumerate(reversed([tuple(range(n_samples))] + S)) :
		#I = np.log(np.take(x, s, axis=axis_smpl))
		#I = [i for i in I if np.isfinite(i)]
		
		
		#if (len(I) < 2) : continue
		#if (min(I) == max(I)) : continue
		
		#I = np.asarray(I)
	
		#t = np.linspace(np.min(I), np.max(I), 100)
		#f = stats.gaussian_kde(I)
		
		#c = ('b' if (n == 0) else 'r')
		#plt.plot(t, f(t), c)
		#if (n != 0) : plt.plot(I, f(I), c + '.')

	#plt.show()
	