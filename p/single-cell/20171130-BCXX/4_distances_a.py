
# RA, 2017-12-07

## ================== IMPORTS :

import re
import math
import pickle
import random
import inspect
import numpy as np

from scipy           import stats
from itertools       import chain
from multiprocessing import cpu_count
from progressbar     import ProgressBar as Progress
from joblib          import Parallel, delayed

## ==================== INPUT :

input_file_BC = "OUTPUT/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

## =================== OUTPUT :

output_file_measurements = "OUTPUT/4_distances_a/measurements.pkl"

## =================== PARAMS :

# Number of parallel computing processes
num_procs = min(math.ceil(cpu_count() / 2), 12)

# Number of dots per graph
dots = 1000

# Number of dimensions to select for each sample (each dot)
DIMS = [10, 100, 1000, 10000, 20000]

metrics = [
	{ 'func' : func, 'data' : data, 'info' : info }
	for data in ['X', 'Z']
	for func in ['eucl', 'corr']
	for info in ['sc', 'fr']
]

metrics.append(
	{ 'func' : 'ks', 'data' : 'X', 'info' : '--' }
)

## ==================== (!!!) :

# Differential expression within a collection of empirical proba
def DE(P) :
	# Distance between two empirical proba
	def dist(p, q) : return stats.ks_2samp(p, q)[0]

	return np.mean([dist(p, q) for p in P for q in P])

# https://en.wikipedia.org/wiki/Silhouette_(clustering)
# D = distance matrix
# S = [[indices of cluster c] for each cluster c]
# Returns the mean silhouette value (score) and
#         the proportion of positive silhouette values (prop)
def SI(D, S) :
	assert(D.shape[0] == D.shape[1])
	def md(i, c) : return np.mean([D[i, j] for j in c])
	A = { c : [md(i, c) for i in c] for c in S }
	B = { c : [min(md(i, d) for d in S if (d != c)) for i in c] for c in S }
	s = { c : [(b - a) / max(b, a) for (a, b) in zip(A[c], B[c])] for c in S }
	#for s in s.values() : print(sorted(s))
	s = list(chain.from_iterable(s.values()))
	return { 'sc' : np.mean(s), 'fr' : sum((i > 0) for i in s) / len(s) }

## ===================== WORK :

# Load the BC data
data = pickle.load(open(input_file_BC, "rb"))
#print(data.keys())

# Expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# Labels of samples of the form BCXX[LN][_Re]_XX 
sample_labels = data['header']

# Remove strange samples
for h in { "BC07LN_20" } :
	X = np.delete(X, sample_labels.index(h), axis=axis_smpl)
	sample_labels.remove(h)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# Get the sample groups
#
# Option 1: BCXX -- by patient
#groups = [n[0:4] for n in sample_labels]
#
# Option 2: BCXX[LN][_Re] -- by batch
groups = [re.findall("(.*)_[0-9]+", n)[0] for n in sample_labels]

# Make groups unique and sort
groups = sorted(list(set(groups)))

# Collect numbers of samples of the form "g"_XX by group
G2S = {
	g : tuple(s for (s, h) in enumerate(sample_labels) if re.match(g + "_[0-9]+", h))
	for g in groups 
}

# Split the expression matrix by groups
G2X = {
	g : np.take(X, G2S[g], axis=axis_smpl) 
	for g in groups 
}

S = list(G2S.values())

Z = stats.mstats.zscore(X, axis=axis_smpl)

def corr(K, data=None, info=None) :
	x = np.take(data, K, axis=axis_gene)
	
	# Covariance & norm products
	C = np.tensordot(x, x, axes=([axis_gene], [axis_gene]))
	V = np.sqrt(np.outer(np.diag(C), np.diag(C)))
	V[V == 0] = 1
	D = 1 - (C / V)
	
	si = SI(D, S)

	return si[info]

def eucl(K, data=None, info=None) :
	x = np.take(data, K, axis=axis_gene)
	
	y = [np.take(x, i, axis=axis_smpl) for i in range(n_samples)]
	D = np.zeros( (n_samples, n_samples) )
	for (i, yi) in enumerate(y) : 
		for (j, yj) in enumerate(y) : 
			D[i, j] = np.linalg.norm(yi - yj)
	
	si = SI(D, S)

	return si[info]

def ks(K) :
	return np.mean([
		DE([np.take(G2X[g], k, axis=axis_gene) for g in groups])
		for k in K
	])

def measure(K, m) :
	if (m['func'] == 'ks') :
		assert(m['data'] == 'X')
		assert(m['info'] == '--')
		return ks(K)
	else :
		func = (corr if (m['func'] == 'corr') else eucl)
		data = (Z if (m['data'] == 'Z') else X)
		info = m['info']
		return func(K, data=data, info=info)

#ks = np.mean([DE([np.take(G2X[g], k, axis=axis_gene) for g in groups]) for k in K])

if (np.inf in DIMS) : DIMS[DIMS.index(np.inf)] = (n_genes - 1)


def main() :
	
	KK = {
		dims : [random.sample(range(n_genes), dims) for _ in range(dots)]
		for dims in DIMS
	}

	# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
	script = inspect.getsource(inspect.getmodule(inspect.currentframe()))

	measurements = { dims : [] for dims in DIMS }
	#
	for dims in DIMS : 
		for metric in metrics :
			print("Aquiring measurements for dims={}, m={}".format(dims, metric))

			measurements[dims].append(
				Parallel(n_jobs=num_procs)(
					delayed(measure)(K, metric)
					for K in Progress()(KK[dims])
				)
			)
			
			# Save results
			pickle.dump(
				{ 
					'measurements' : measurements, 
					'DIMS'         : DIMS, 
					'metrics'      : metrics,
					'script'       : script
				}, 
				open(output_file_measurements, "wb")
			)

if (__name__ == "__main__") :
	main()

