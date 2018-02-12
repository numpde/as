
# RA, 2018-02-12

## ================== IMPORTS :

import numpy as np
import pickle
import math
import os.path
import inspect
import sys
import pandas as pd

from collections import defaultdict
from scipy import stats
from itertools   import chain
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	# Extract the list of relevant genes from here
	'BC data' : "OUTPUT/0_prepare/UV/bcxx.pkl",
	
	# 
	'GO=>Symb' : "ORIGINALS/go-1/go2symb.txt",
	'GO=>Info' : "ORIGINALS/go-1/go2name.txt",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	# 
	'GO=>CI' : "OUTPUT/0_go2ci/UV/go2ci.pkl",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## =================== PARAMS :

# Test mode
TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' : TESTMODE or min(12, math.ceil(cpu_count() / 1.5)),
	
	# Category size window width for the windowed quantiles
	'window width' : 2,
}

## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
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


## ===================== WORK :



#[ LOAD BC DATASET ]#

# Load the BC data
BCXX = pickle.load(open(IFILE['BC data'], 'rb'))

# Expression matrix
X = BCXX['X']

# Data layout: genes x samples
(axis_gene, axis_smpl) = (0, 1)

# Drop unexpressed genes
X = X[ X.sum(axis=axis_smpl) != 0 ]

# H2I : HGNC Symbol --> Gene indices in BC data
H2I = dict(zip(X.index, range(len(X.index))))

# Clusters/groups
G2S = { 
	g : tuple(s for (s, h) in SH)
	for (g, SH) in BCXX['B2SH'].items() 
}
S = list(G2S.values())

#[ LOAD GO TERMS & ANNOTATION ]#

# GO2H : GO ID --> [Gene HGNC symbols]
GO2H = {
	r[0] : r[1].split('|')
	for r in pd.read_csv(IFILE['GO=>Symb'], sep='\t', index_col=0).itertuples()
}
# Filter out those symbols that not in the expression table
GO2H = { go : [h for h in H if (h in H2I)] for (go, H) in GO2H.items()}

# GO2T : GO ID --> GO category name
GO2T = {
	r[0] : r[1]
	for r in pd.read_csv(IFILE['GO=>Info'], sep='\t', index_col=0).itertuples()
}


# All GO terms
GO = sorted(GO2T.keys())

if TESTMODE : GO = sorted(np.random.permutation(GO)[0:100])


#[ GO TERMS vs BC DATASET ]#

#for go in GO2T.keys() :
	#if (go not in GO2H) :
		##print("NOTE: {} has no gene list".format(go))
		#GO2H[go] = set()

# GO2I : GO ID --> Gene indices in BC data
GO2I = {
	go : np.asarray([H2I[h] for h in H])
	for (go, H) in GO2H.items()
}


#[ PREPARE DATA ]#

Z = np.moveaxis(stats.mstats.zscore(X.as_matrix(), axis=axis_smpl), axis_gene, 0)

#[ COMPUTE THE CLUSTERING INDEX FOR EACH GO TERM ]#

def goci_from_go(go) :
	if (0 == len(GO2I.get(go, []))) : return (go, None)
	
	return ( go, CI(Z[GO2I[go], :], 0) )

print("Computing clustering indices")

GO2CI = dict(
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(goci_from_go)(go)
		for go in Progress()(np.random.permutation(GO)) 
	)
)

# N2CI : size of GO term --> [clustering indices]
N2CI = defaultdict(list)
for (go, I) in GO2I.items() : 
	if (go in GO2CI) and len(I) :
		N2CI[len(I)].append(GO2CI[go])


#[ COMPUTE THE WINDOWED QUANTILES FOR EACH GO TERM ]#

def gowq_from_go(go) :
	I = GO2I.get(go, [])
	
	if (len(I) == 0) : return (go, None)

	w = math.sqrt(PARAM['window width'])

	p = stats.percentileofscore(
		[ci for (n, CI) in N2CI.items() for ci in CI if (len(I)/w <= n <= len(I)*w)],
		GO2CI[go]
	)
	
	q = min(1, max(0, p / 100))
	
	return (go, q)

print("Computing windowed quantiles")

GO2WQ = dict(
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(gowq_from_go)(go)
		for go in Progress()(np.random.permutation(GO)) 
	)
)


#[ SAVE ]#

pickle.dump(
	{ 
		'GO2H'   : GO2H,
		'GO2T'   : GO2T,
		'GO2CI'  : GO2CI,
		'N2CI'   : N2CI,
		'GO2WQ'  : GO2WQ,
		'G2S'    : G2S,
		'S'      : S, 
		'PARAM'  : PARAM,
		'script' : THIS,
	}, 
	open(OFILE['GO=>CI'], "wb")
)
