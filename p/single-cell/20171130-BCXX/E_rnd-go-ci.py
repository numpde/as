
# RA, 2018-01-16

# Run as
#	python3 E_rnd-go-ci.py COMPUTE
#	python3 E_rnd-go-ci.py LIST
#	python3 E_rnd-go-ci.py PLOT

## ================== IMPORTS :

import re
import os
import sys
import math
import time
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
from glob            import glob as list_files
from fractions       import Fraction

from sklearn.metrics.pairwise import euclidean_distances as euc_dist

# https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
from scipy.stats     import gaussian_kde

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_prepare/UV/bcxx.pkl",
	'GO=>CI'   : "OUTPUT/0_go2ci/UV/go2ci.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'runs' : "OUTPUT/E_rnd-go-ci/UV/E_rnd-go-ci_{ver}.pkl",
	
	'list' : "OUTPUT/E_rnd-go-ci/list_pivot={pivot}_gonly={union}.txt",
	
	'plot' : "OUTPUT/E_rnd-go-ci/freq_pivot={pivot}.pdf",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' : TESTMODE or min(12, math.ceil(cpu_count() / 1.2)),
	
	# Number of random subsets per GO category
	'M' : 10,
	
	# Exclude genes w/o GO category
	'GO union only' : True,
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
BCXX = pickle.load(open(IFILE['BC data'], 'rb'))

# Expression matrix
BC_X = BCXX['X']

# Data layout: genes x samples
(axis_gene, axis_smpl) = (0, 1)

## Number of samples / genes in the expression matrix
#(n_samples, n_genes) = (BC_X.shape[axis_smpl], BC_X.shape[axis_gene])

# HGNC Symbols in the table
BC_H = list(BC_X.index)

# H2I : HGNC Symbol --> Gene indices in BC data
H2I = dict(zip(BC_X.index, range(len(BC_X.index))))

# Clusters/groups
G2S = { 
	g : tuple(s for (s, h) in SH)
	for (g, SH) in BCXX['B2SH'].items() 
}
S = list(G2S.values())

#[ GO DATA ]#

# Clustering indices data bundle
CI_data = pickle.load(open(IFILE['GO=>CI'], 'rb'))

# GO2CI : GO ID --> Clustering index
GO2CI = CI_data['GO2CI']

# GO2H : GO ID --> [HGNC Symbols in the table]
GO2H = CI_data['GO2H']

# GO2T : GO ID --> GO category name
GO2T = CI_data['GO2T']

# GO2WQ : GO ID --> clustering index windowed quantile
GO2WQ = CI_data['GO2WQ']


## =============== PREPROCESS :

#[ Remove repeated GO categories ]#

H2GO = defaultdict(set)
for (go, H) in GO2H.items() : H2GO[hash('+'.join(sorted(H)))].add(go)
H2GO = dict(H2GO)
#
# Check for no hash collisions
assert(all((1 == len(set('+'.join(sorted(GO2H[go])) for go in GO))) for GO in H2GO.values()))
#
# Non-redundant GO categories and their aliases
GO2A = { min(GO) : sorted(GO) for GO in H2GO.values() }
#
#print("{} of {} GO categories are non-redundant".format(len(GO2A), len(GO2H)))
assert(10000 <= len(GO2A) <= 30000), "Unexpected number of GO terms"
#
del H2GO

if TESTMODE : 
	PARAM['M'] = 2
	GO2A = dict(list(GO2A.items())[0:200])

def restrict(GO2X, GO) :
	return { go : x for (go, x) in GO2X.items() if (go in GO) }

GO2H  = restrict(GO2H, GO2A.keys())
GO2CI = restrict(GO2CI, GO2A.keys())
GO2WQ = restrict(GO2WQ, GO2A.keys())


## =========== COMPUTING WORK :

# Exclude genes w/o GO category
if PARAM['GO union only'] :
	H = set(chain.from_iterable(GO2H.values()))
	BC_X = np.take(BC_X, [n for (n, h) in enumerate(BC_H) if h in H], axis=axis_gene)
	print("Keeping {}/{} genes".format(len(H), len(BC_H)))
	# Those are invalidated:
	del H2I
	del BC_H

Z = stats.mstats.zscore(BC_X.as_matrix(), axis=axis_smpl)
assert(axis_gene == 0) # For slicing the Z matrix

def job(I) : return ( len(I), CI(Z[I, :], axis_gene) )

def COMPUTE() :
	
	# For each GO category, generate M random subsets of the same size
	M = PARAM['M']
	
	ver = time.strftime("%Y%m%d-%H%M%S")
	print("Results will be written to:")
	print(OFILE['runs'].format(ver=ver))
	
	print("Preparing random subsets")
	II = list(
		random.sample(range(len(Z)), len(H))
		for (_, H) in Progress()(GO2H.items())
		for _ in range(M)
		if len(H)
	)
	
	# Main computation loops
	print("Computing clustering indices")
	NC = list(
		Parallel(n_jobs=PARAM['#proc'])(
			delayed(job)(I)
			for I in Progress()(II)
		)
	)
	
	# Save the results
	pickle.dump(
		{ 'NC' : NC, 'script' : THIS, 'PARAM' : PARAM },
		open(OFILE['runs'].format(ver=ver), 'wb')
	)


## ======== LIST GO QUANTILES :

def LIST() :
	
	# Load the clustering indices of random subsets
	RUNS = [
		pickle.load(open(f, 'rb'))
		for f in sorted(list_files(OFILE['runs'].format(ver="*")))
	]
	
	# Fix old format (before 2018-01-23)
	for (i, run) in enumerate(RUNS) :
		if not ('PARAM' in run) :
			RUNS[i]['PARAM'] = { 'GO union only' : False }
	
	# Load the clustering indices of random subsets
	NC = list(chain.from_iterable(
		run['NC']
		for run in RUNS
		if (run['PARAM']['GO union only'] == PARAM['GO union only'])
	))
	
	Q = dict()
	for K in [2**k for k in range(0, 13)] :
		print("Pivot:", K)
		
		w = math.sqrt(2)
		C = [c for (n, c) in NC if (K/w <= n < K*w)]
		
		QK = []
		
		for (go, H) in Progress()(GO2H.items()) :
			if not (K/w <= len(H) < K*w) : continue
			if not GO2CI.get(go, None) : continue
		
			ci = GO2CI[go]
			
			q = np.mean([(c < ci) for c in C])
			
			QK.append( (go, len(H), K, GO2T[go], GO2WQ[go], q, len(C)) )
		
		Q[K] = sorted(QK, key=(lambda x : x[5]))
	
	filename = OFILE['list'].format(pivot="all", union=PARAM['GO union only'])
	
	with open(filename, 'w') as f :
		print("GO term", "GO size", "Size pivot", "GO name", "CI quantile (GO)", "CI quantile (random)", "Number of random subsets", sep='\t', file=f)
		for (K, QK) in sorted(Q.items()) :
			for q in QK : 
				print(*q, sep='\t', file=f)


## ============ PLOTTING WORK :

def PLOT() :
	
	# Load the clustering indices of random subsets
	RUNS = [
		pickle.load(open(f, 'rb'))
		for f in sorted(list_files(OFILE['runs'].format(ver="*")))
	]
	
	# Fix old format (before 2018-01-23)
	for (i, run) in enumerate(RUNS) :
		if not ('PARAM' in run) :
			RUNS[i]['PARAM'] = { 'GO union only' : False }
	
	# Add the tag according to underlying sampling set
	for (i, run) in enumerate(RUNS) :
		RUNS[i]['tag'] = ('m' if run['PARAM']['GO union only'] else 'r')
	
	# Tag to legend
	T2L = {
		'r' : "Random subsets (all genes)", 
		'm' : "Random subsets (GO genes)", 
		'g' : "GO categories",
		'-' : "5% and 95%",
	}
		
	# RUNS[i]['NC'] is a list of the form
	#    [(n1, c1), (n2, c2), ...]
	# where n is the subset size and c is the clustering index
	
	# Make a similar list for the GO categories
	NC_GO = [(len(E), GO2CI[go]) for (go, E) in GO2H.items() if (len(E) and (go in GO2CI))]
	
	# Collect all those lists with a tag
	#    'r' = random subset
	#    'm' = random subset (genes from the GO union)
	#    'g' = GO category
	NC_ALL = [(runs['NC'], runs['tag']) for runs in RUNS] + [(NC_GO, 'g')]
	
	for K in [2**k for k in range(0, 11)] :
		print("Subset size:", K)
		
		w = math.sqrt(2)

		plt.clf()
		
		# Handles and legends will be filled in here
		H = { tag : [] for tag in T2L.keys() }
		L = { tag : [] for tag in T2L.keys() }
		
		for (NC, tag) in NC_ALL :
			
			nc = [(n, c) for (n, c) in NC if (K/w <= n < K*w)]
			
			if not nc : continue
			
			# Filter subsets by size
			(N, C) = zip(*nc)
			
			f = gaussian_kde(C)
			
			t = np.linspace(-1, 1, 257)
			
			plt.plot(t, f(t), ('-' + tag))
			
			h = plt.plot(-1, 0, ('-' + tag), linewidth=4)[0]
			
			H[tag].append( h )
			L[tag].append( "{} of size {}--{}".format(T2L[tag], min(N), max(N)) )
			
			# Indicate the 5% and 95% sections
			for p in np.percentile(C, [5, 95]) :
				plt.plot((p, p), (0, f(p)), ('--' + tag))
		
			# Plot also the data points for the GO categories
			if (tag == 'g') :
				plt.plot(C, f(C), ('.' + tag), markersize=3)
		
		H['-'] = [ plt.plot(-1, 0, '--k')[0] ]
		L['-'] = [ T2L['-'] ]
		
		ylim = max(plt.ylim())
		plt.yticks(range(15))
		
		# Optional:
		ylim = 12.5
		
		xx = np.linspace(-1, 1, 9)
		plt.xticks(xx, [Fraction(x) for x in xx])
		
		plt.xlim([-1, 1])
		plt.ylim([0, ylim])
		
		plt.legend(
			[ H[tag][0] for tag in ['r', 'm', 'g', '-'] if (tag in H) ],
			[ L[tag][0] for tag in ['r', 'm', 'g', '-'] if (tag in L) ],
			#loc = ('upper left' if (np.median(C) >= 0) else 'upper right'),
			loc = 'upper left'
		)
		
		plt.xlabel("Clustering index")
		plt.ylabel("Relative frequency")
		
		plt.tight_layout()
		
		plt.savefig(OFILE['plot'].format(pivot=K))


## ===================== MAIN :

if (__name__ == "__main__") :
	
	if ("COMPUTE" in sys.argv) : COMPUTE()
	if ("LIST"    in sys.argv) : LIST()
	if ("PLOT"    in sys.argv) : PLOT()

