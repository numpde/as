
# RA, 2017-12-20

# 

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import numpy as np

from scipy           import stats
from itertools       import chain
from collections     import Counter
from progressbar     import ProgressBar as Progress

from multiprocessing import cpu_count
from joblib          import Parallel, delayed

## ==================== INPUT :

IFILE = {
	"BC data" : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	"GO -> ENSG" : "OUTPUT/0_e2go/go2e.txt",
}

## =================== OUTPUT :

OFILE = {
	"Results" : "OUTPUT/6_low-high-end_b/UV/results.pkl",
}

## =================== PARAMS :

# Test mode
TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of random subsets
	'#dots' : (1000 if not TESTMODE else 22),
	
	# Number of dimensions per random subset
	'#dims' : [5, 10, 15],
	
	# GO ID should appear in that many genes to be considered
	'e/go cut-off' : 20,
	
	# Number of parallel computing processes
	'#proc' : min(24, math.ceil(cpu_count() / 2)),
}

## ==================== PREPA :

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)
	
# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

## ===================== WORK :

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

def job(X, axis, dims, S) :
	# Random subset of genes of size dims
	K = random.sample(range(X.shape[axis]), dims)
	# Clustering index
	i = silhouette(cos_dist(np.take(X, K, axis=axis), axis), S)
	i = list(chain.from_iterable(i.values()))
	if i :
		return np.mean(i) # OR: np.mean((x > 0) for x in i)
	else :
		return None

def main() :

	# [ LOAD GO TERMS ]

	# GO2E : GO ID --> [ENSG IDs]
	GO2E = dict(
		(go, E)
		for (go, E) in [
			(go_E[0], go_E[1:])
			for go_E in [
				L.rstrip().split('\t') 
				for L in open(IFILE["GO -> ENSG"], 'r')
			]
		]
		if (len(E) >= PARAM['e/go cut-off'])
	)
	
	if TESTMODE : GO2E = dict(list(GO2E.items())[0:7])
	
	print("Profiling {} GO IDs".format(len(GO2E)))
	
	# [ LOAD BC DATASET ]
	
	# Load the BC data
	BC_data = pickle.load(open(IFILE["BC data"], "rb"))

	# Expression matrix
	X = BC_data['X']
	
	# ENSG IDs
	E = BC_data['gene_id']

	# Labels for axis/dimension of BC data
	(axis_smpl, axis_gene) = (BC_data['axis_smpl'], BC_data['axis_gene'])

	# Z transform
	Z = stats.mstats.zscore(X, axis=axis_smpl)
	
	# Prevent accidental use of the original:
	del X 

	# Labels of samples of the form BCXX[LN][_Re]_XX 
	sample_labels = BC_data['header']
	
	# Number of samples / genes in the expression matrix
	(n_samples, n_genes) = (Z.shape[axis_smpl], Z.shape[axis_gene])
	
	# S = [(indices of items in c) for each cluster c]
	S = list(tuple(s for (s, _) in SH) for SH in BC_data['B2SH'].values())
	
	# [ RANDOM SUBSETS ]
	
	GO2I = dict()
	for dims in PARAM['#dims'] :
		
		# GO -> [clustering indices]
		GO2I[dims] = dict( (go, []) for go in GO2E.keys() )
		
		for (go, E) in GO2E.items() :
			
			print("Profiling '{}' by subsets of size {}".format(go, dims))
			
			# Iterate over random subsets K
			
			GO2I[dims][go] = Parallel(n_jobs=PARAM['#proc'])(
				delayed(job)(Z, axis_gene, dims, S)
				for _ in Progress()(range(PARAM['#dots']))
			)
			
			# Do not save results if in test mode
			if TESTMODE : continue
		
			pickle.dump(
				{'GO2I' : GO2I, 'script' : THIS}, 
				open(OFILE["Results"], "wb")
			)

if (__name__ == "__main__") :
	main()
