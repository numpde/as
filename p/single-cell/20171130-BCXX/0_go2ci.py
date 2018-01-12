
# RA, 2018-01-04

## ================== IMPORTS :

import numpy as np
import pickle
import math
import os.path
import inspect
import sys

from collections import defaultdict
from scipy import stats
from itertools   import chain
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	# Extract the list of relevant genes from here
	'BC data' : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	
	# 
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	
	#
	'GO=>Info' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
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

PARAM = {
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.5)),
	
	# Window width for the windowed quantiles
	'window width' : 2,
}

# Test mode
TESTMODE = ("TEST" in sys.argv)

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
def cos_dist(X, axis) :
	# Covariance & norm products
	C = np.tensordot(X, X, axes=([axis], [axis]))
	V = np.sqrt(np.outer(np.diag(C), np.diag(C)))
	V[V == 0] = 1
	D = 1 - (C / V)
	return D

## ===================== WORK :



#[ LOAD BC DATASET ]#

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

## E2X : BC ENSG --> BC X data
#E2X = dict(zip(BC_E, np.moveaxis(X, axis_gene, 0)))

# E2I : BC ENSG --> Gene indices in BC data
E2I = dict(zip(BC_E, range(len(BC_E))))

# Clusters/groups
G2S = { 
	g : tuple(s for (s, h) in SH)
	for (g, SH) in BC_data['B2SH'].items() 
}
S = list(G2S.values())

#[ LOAD GO TERMS & ANNOTATION ]#

# GO2E : GO ID --> [ENSG IDs]
GO2E = {
	go_E[0] : set(go_E[1:])
	for go_E in [
		L.rstrip().split('\t') 
		for L in open(IFILE['GO=>ENSG'], 'r')
	]
}

# GO2T : GO ID --> GO category name
GO2T = {
	go : data['name']
	for (go, data) in pickle.load(open(IFILE['GO=>Info'], 'rb')).nodes(data=True)
}

# All GO terms
GO = list(GO2T.keys())

if TESTMODE : GO = GO[0:100]


#[ GO TERMS vs BC DATASET ]#

for go in GO2T.keys() :
	if (go not in GO2E) :
		#print("NOTE: {} has no gene list".format(go))
		GO2E[go] = set()

# GO2I : GO ID --> Gene indices in BC data
GO2I = {
	go : np.asarray([E2I[e] for e in E])
	for (go, E) in GO2E.items()
}


#[ PREPARE DATA ]#

Z = np.moveaxis(stats.mstats.zscore(X, axis=axis_smpl), axis_gene, 0)
del X

assert(Z.shape[0] == n_genes)
assert(Z.shape[1] == n_samples)


#[ COMPUTE THE CLUSTERING INDEX FOR EACH GO TERM ]#

def goci_from_go(go) :
	if (len(GO2I[go]) == 0) : return (go, None)

	return (
		go, 
		np.mean([
			np.sign(x) 
			for x in chain.from_iterable(silhouette(cos_dist(Z[GO2I[go], :], 0), S).values())
		])
	)

print("Computing clustering indices")

GO2CI = dict(
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(goci_from_go)(go)
		for go in Progress()(np.random.permutation(GO)) 
	)
)

# N2CI : size of GO term --> [clustering indices]
N2CI = defaultdict(list)
for (go, E) in GO2E.items() : 
	if (go in GO2CI) and len(E) :
		N2CI[len(E)].append(GO2CI[go])


#[ COMPUTE THE WINDOWED QUANTILES FOR EACH GO TERM ]#

def gowq_from_go(go) :
	E = GO2E[go]
	
	if (len(E) == 0) : return (go, None)

	w = math.sqrt(PARAM['window width'])

	p = stats.percentileofscore(
		[ci for (n, CI) in N2CI.items() for ci in CI if (len(E)/w <= n <= len(E)*w)],
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
		'GO2E'   : GO2E,
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
