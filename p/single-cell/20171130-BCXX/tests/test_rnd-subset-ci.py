
# RA, 2018-01-16

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.sparse    import lil_matrix
from collections     import defaultdict
from string          import ascii_lowercase
from numpy.matlib    import repmat
from scipy           import stats
from scipy.constants import golden as phi
from itertools       import chain
from multiprocessing import cpu_count
from joblib          import Parallel, delayed
from progressbar     import ProgressBar as Progress

from sklearn.metrics.pairwise import euclidean_distances as euc_dist

## ==================== INPUT :

IFILE = {
	'BC data'  : "../OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO graph' : "../OUTPUT/0_go-graph/UV/go-graph.pkl",
	'GO=>CI'   : "../OUTPUT/0_go2ci/UV/go2ci.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'runs' : "./test_rnd-subset-ci.pkl",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	# Figure formats
	'ext' : ['png', 'pdf'],
	
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.2)),
}

mpl.rcParams['axes.labelsize'] = 'large'

## ====================== AUX :

# This script
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

## ====================== (!) :


# https://en.wikipedia.org/wiki/Silhouette_(clustering)
# D = distance matrix
# S = [[indices of cluster c] for each cluster c]
# Returns the silhouette values by cluster
def silhouette(D, S) :
	assert(D.shape[0] == D.shape[1])
	def md(i, c) : return np.mean([D[i, j] for j in c])
	A = { c : [md(i, c) for i in c] for c in S }
	B = { c : [min(md(i, d) for d in S if (d != c)) for i in c] for c in S }
	s = { c : [(b - a) / max(b, a) for (a, b) in zip(A[c], B[c]) if max(b, a)] for c in S }
	#for s in s.values() : print(sorted(s))
	return s

# Compute a distance matrix as (1 - cos(angle))
def cos_dist(X) :
	# Covariance & norm products
	C = np.tensordot(X, X, axes=([1], [1]))
	V = np.sqrt(np.outer(np.diag(C), np.diag(C)))
	V[V == 0] = 1
	D = 1 - (C / V)
	return D

# Clustering index
def CI(X, feature_axis) :

	def dist(X) : return cos_dist(np.moveaxis(X, feature_axis, 1))
	#def dist(X) : return euc_dist(np.moveaxis(X, feature_axis, 1))

	return np.mean([
		np.sign(x) 
		for x in chain.from_iterable(silhouette(dist(X), S).values())
	])

## ===================== DATA :

#[ BC DATA ]#

# Load the BC data
BC_data = pickle.load(open(IFILE['BC data'], 'rb'))

# Expression matrix
BC_X = BC_data['X']

# Rearrange data axes
(axis_gene, axis_smpl) = (0, 1)
BC_X = np.moveaxis(BC_X, BC_data['axis_gene'], axis_gene)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (BC_X.shape[axis_smpl], BC_X.shape[axis_gene])

# ENSG IDs
BC_E = BC_data['gene_id']
assert(len(BC_E) == BC_X.shape[axis_gene]), "Inconsistent gene info"

# E2I : BC ENSG --> Gene indices in BC data
E2I = dict(zip(BC_E, range(len(BC_E))))

# Clusters/groups
G2S = { 
	g : tuple(s for (s, h) in SH)
	for (g, SH) in BC_data['B2SH'].items() 
}
S = sorted(G2S.values())

#[ GO DATA ]#

# Clustering indices data bundle
CI_data = pickle.load(open(IFILE['GO=>CI'], 'rb'))

# GO2E : GO ID --> Clustering index
GO2CI = CI_data['GO2CI']

# GO2E : GO ID --> [ENSG IDs]
GO2E = CI_data['GO2E']

# GO2I : GO ID --> Gene indices in BC data
GO2I = {
	go : np.asarray([E2I[e] for e in E])
	for (go, E) in GO2E.items()
}

# GO2T : GO ID --> GO category name
GO2T = CI_data['GO2T']

# GO2WQ : GO ID --> clustering index windowed quantile
GO2WQ = CI_data['GO2WQ']


## =============== PREPROCESS :

#[ Remove repeated GO categories ]#

H2GO = defaultdict(set)
for (go, E) in GO2E.items() : H2GO[hash('+'.join(sorted(E)))].add(go)
H2GO = dict(H2GO)
#
# Check for no hash collisions
assert(all((1 == len(set('+'.join(sorted(GO2E[go])) for go in GO))) for GO in H2GO.values()))
#
# Non-redundant GO categories and their aliases
GO2A = { min(GO) : sorted(GO) for GO in H2GO.values() }
#
#print("{} of {} GO categories are non-redundant".format(len(GO2A), len(GO2E)))
assert(10000 <= len(GO2A) <= 30000), "Unexpected number of GO terms"
#
del H2GO

def restrict(GO2X, GO) :
	return { go : x for (go, x) in GO2X.items() if (go in GO) }

GO2E  = restrict(GO2E, GO2A.keys())
GO2I  = restrict(GO2I, GO2A.keys())
GO2CI = restrict(GO2CI, GO2A.keys())
GO2WQ = restrict(GO2WQ, GO2A.keys())


## ===================== WORK :

Z = stats.mstats.zscore(BC_X, axis=axis_smpl).tolist()

def job(N) :
	return CI(np.vstack(random.sample(Z, N)), axis_gene) 

def compute() :
	assert(axis_gene == 0)
	
	M = 1000 # Number of subsets per size
	NN = []
	CC = []
	for n in range(0, 20) :
		N = 2 ** n
		if (N > n_genes) : break
	
		print("Subset size:", N)
		
		C = set(
			Parallel(n_jobs=PARAM['#proc'])(
				delayed(job)(N)
				for _ in Progress()(range(M))
			)
		)
		
		NN.append(N)
		CC.append(C)
		
		pickle.dump(
			{ 'NN' : NN, 'CC' : CC, 'M' : M },
			open(OFILE['runs'], 'wb')
		)
		
		#plt.ion()
		#plt.semilogx(N, C, '.k', markersize=1)
		#plt.pause(0.1)
		#plt.show()
	
	#input()

###

if (__name__ == "__main__") :
	compute()

