
# RA, 2017-12-12

# Empirical frequency of GO terms by clustering index
# over random samples of genes.

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import numpy as np
#import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
from scipy.stats     import gaussian_kde

from scipy           import stats
from itertools       import chain
from collections     import Counter
from progressbar     import ProgressBar as Progress

#from multiprocessing import cpu_count
# Number of parallel computing processes
#'#proc' : math.ceil(cpu_count() / 2),

## ==================== INPUT :

IFILE = {
	"BC data" : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	"ENSG to GO" : "OUTPUT/0_e2go/e2go.txt"
}

## =================== OUTPUT :

OFILE = {
	"Results" : "OUTPUT/6_low-high-end_a/UV/results.pkl",
	#"Frequency" : "OUTPUT/6_low-high-end_a/go-freq_dims={dims}.{ext}"
}

## =================== PARAMS :

PARAM = {
	# Number of random subsets
	'#dots' : (100000 if not ("TEST" in sys.argv) else 22),
	
	# Number of dimensions per random subset
	'#dims' : [10, 100, 1000]
}

## ==================== PREPA :

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)
	
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
	s = { c : [(b - a) / max(b, a) for (a, b) in zip(A[c], B[c])] for c in S }
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

def main() :

	# [ LOAD GO TERMS ]

	# Read ENSG-GO associations from file
	# E2GO : ENSG ID --> [GO IDs]
	with open(IFILE["ENSG to GO"], 'r') as f :
		E_GO = [L.rstrip().split('\t') for L in f]
		E2GO = { e_go[0] : e_go[1:] for e_go in E_GO }
		del E_GO
	
	# All known GO IDs
	GO = list(set(chain.from_iterable(E2GO.values())))
	
	# GO2E : GO ID --> [ENSG IDs]
	GO2E = { go : [] for go in GO }
	for (e, gos) in E2GO.items() : 
		for go in gos :
			GO2E[go].append(e)

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
	
	GO2I = {}
	for dims in PARAM['#dims'] :
		
		print("Computing with gene subsets of size {}".format(dims))
		
		# GO to index
		GO2I[dims] = { go : [] for go in GO }
		
		# Iterate over random subsets K
		for _ in Progress()(range(PARAM['#dots'])) :
			
			# Random subset of genes of size dims
			K = random.sample(range(n_genes), dims)
			
			# Clustering index
			i = silhouette(cos_dist(np.take(Z, K, axis=axis_gene), axis_gene), S)
			i = list(chain.from_iterable(i.values()))
			i = np.mean(i) # OR: i = np.mean((x > 0) for x in i)
			
			Ek = [E[k] for k in K if (E[k] in E2GO)]
			
			go2c = Counter(go for e in Ek for go in E2GO[e])
			
			for (go, c) in go2c.items() :
				# Expected number of occurrences of this GO ID in the sample
				x = len(GO2E[go]) / n_genes * dims
				
				if (c >= x) : GO2I[dims][go].extend([i] * c)
		
		# Absolute frequency of the most frequent GO ID
		maxI = max(len(I) for I in GO2I[dims].values())
		
		pickle.dump({'GO2I' : GO2I}, open(OFILE["Results"], "wb"))
		
		#for (go, I) in GO2I[dims].items() :
			#if (len(I) == 0) : continue
			#if (len(I) < 0.1 * maxI) : continue
			#if (min(I) == max(I)) : continue

			#t = np.linspace(min(I), max(I), 100)
			#f = gaussian_kde(I)(t)
			
			#plt.plot(t, f)
			
			##plt.hist(I, histtype='step', density=True, stacked=True)
		
		#plt.xlabel("Clustering index")
		#plt.ylabel("Relative empirical frequency")
		#plt.savefig(OFILE["Frequency"].format(dims=dims, ext="png"))

if (__name__ == "__main__") :
	main()

