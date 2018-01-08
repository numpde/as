
# RA, 2017-12-20

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import networkx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy import stats
from scipy.constants import golden as phi
from itertools import chain
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'GO=>Info' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'gene-ks'  : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'2d-plot' : "OUTPUT/8_de-vs-ci/de-vs-ex_{go}_{log}.{ext}",
	'2d-info'  : "OUTPUT/8_de-vs-ci/de-vs-ex_{go}_info.txt",
	
	'ci-vs-de' : "OUTPUT/8_de-vs-ci/ci-vs-de_{log}_{way}.{ext}",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	# GO terms of interest
	'GO filter' : {
		"GO:0001525", # angiogenesis
		"GO:0006281", # DNA-repair
		"GO:0006955", # immune response
		"GO:0007049", # cell cycle
		"GO:0016477", # cell migration
		"GO:0004984", # olfactory receptor activity
	},
	
	# Figure formats
	'ext' : ['png'],
	
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.2)),
}

mpl.rcParams['axes.labelsize'] = 'large'

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

# E2X : BC ENSG --> BC X data
E2X = dict(zip(BC_E, np.moveaxis(X, axis_gene, 0)))


#[ LOAD KS DIFFERENTIAL EXPRESSION DATA ]#

KS_data = pickle.load(open(IFILE['gene-ks'], 'rb'))
#print(KS.keys())

KS_OP = 1
KS_meta = KS_data['KS_meta'][KS_OP]
assert(KS_meta in ['max', 'mean', 'median', 'min', 'std'])

# Differential expression by ENSG
E2DE = { e : ks[KS_OP] for (e, ks) in KS_data['E2KS'].items() }

# Group to [sample indices]
G2S = KS_data['G2S']

# Clusters/groups by group
G = sorted(list(G2S.keys()))
S = [ tuple(G2S[g]) for g in G ]
del G2S # Use G and S


#[ COMPUTE EXPRESSION MAGNITUDE ]#

E2EX = { e : np.mean(E2X[e]) for e in BC_E }


#[ LOAD GO TERMS & ANNOTATION ]#

# GO2E : GO ID --> [ENSG IDs]
# Read it from file
GO2E = {
	go_E[0] : set(go_E[1:])
	for go_E in [
		L.rstrip().split('\t') 
		for L in open(IFILE['GO=>ENSG'], 'r')
	]
}

# GO2T : GO ID --> GO category name
#
GO2T = {
	n : data['name']
	for (n, data) in pickle.load(open(IFILE['GO=>Info'], 'rb')).nodes(data=True)
}
#
# Old version, reading from text file
#
#GO2T = dict(
	#tuple(L.rstrip().split('\t')[:2])
	#for L in open(IFILE['GO=>Info'], 'r').readlines()[1:]
#)

for go in list(PARAM['GO filter']) :
	if (go is not None) and (go not in GO2E) :
		print("NOTE: {} has no gene list".format(go))
		PARAM['GO filter'].remove(go)

# Narrow down the gene selection (to each GO term of interest)
GO2E = {
	go : [
		BC_E.index(e) 
		for e in sorted(set(BC_E) & GO2E.get(go), key=(lambda e : -E2DE[e]))
	]
	for go in sorted(PARAM['GO filter'])
}
	
#[ PLOT CI VS DE-CUTOFF ]#

if True :
	
	np.random.seed(0)
	
	def clustering_score(X) :
		return np.mean([np.sign(x) for x in chain.from_iterable(silhouette(cos_dist(X, 0), S).values())])
	
	xlabel = { 
		'm' : "Number of genes in order of mean expression",
		'+' : "Number of most DE genes",
		'-' : "Number of least DE genes",
		'r' : "Number of genes in a random order",
	}

	for way in ['m', 'r', '+', '-'] :
		
		plt.figure(figsize=(math.ceil(6*phi), 6), dpi=150)
		
		GOXS = []
		
		for (go, E) in GO2E.items() :
			print(way, "/", go)
			
			Y = stats.mstats.zscore(np.moveaxis(X, axis_gene, 0)[E], axis=1)
			V = Y[np.random.permutation(len(Y))]
			W = sorted(Y, key=(lambda y : -np.mean(y)))
			
			def select(n, way) :
				if (way == 'r') : return V[0:n]
				if (way == 'm') : return W[0:n]
				if (way == '+') : return Y[0:n]
				if (way == '-') : return Y[-n:]
			
			# Compute silhouette values
			s = Parallel(n_jobs=PARAM['#proc'])(
				delayed(clustering_score)(select(n, way))
				for n in Progress()(range(1, len(E)))
			)
			
			#x = np.linspace(0.0, 1.0, len(E))[1:].tolist() # relative
			x = list(range(1, len(E))) # absolute
			
			GOXS.append((go, x, s))
		
		L = [] # Legend
		for (go, x, s) in sorted(GOXS, key=(lambda x : -np.max(x[2]))) :
			plt.plot(x, s)
			L.append("{} ({})".format(GO2T[go], len(GO2E[go])))
		
		plt.ylim(-1, +1)
		
		plt.xlabel(xlabel[way])
		plt.ylabel("Clustering index")
		
		plt.legend(L, loc='upper left')
		
		for scale in ['linear', 'log'] :
			plt.xscale(scale)
			for ext in PARAM['ext'] :
				plt.savefig(OFILE['ci-vs-de'].format(log=scale[0:3], way=way, ext=ext))
		
		plt.close()

if ("ONLY1" in sys.argv) : exit()

#[ FOR EACH GO TERM OF INTEREST ... ]#

for (go, E) in GO2E.items() :
	
	# GO ID for file names
	go_safe = str(go).replace(':', '-')
	
	#[ ... PLOT 2-GENE EXPRESSION ]#
	
	# Pick the top 2 for a 2d plot
	E = E[0:2]
	
	plt.figure(figsize=(math.ceil(6*phi), 6), dpi=150)
	plt.xscale('linear')
	plt.yscale('linear')
	
	colors = plt.get_cmap('hsv')(np.linspace(0.1, 1.0, len(G))).tolist()
	
	L = [] # Legend
	for (g, s) in zip(G, S) :
		Y = np.take(X, E, axis=axis_gene)
		Y = np.take(Y, s, axis=axis_smpl)
		a = np.take(Y, 0, axis=axis_gene)
		b = np.take(Y, 1, axis=axis_gene)
		plt.plot(a, b, '.', markersize=8, color=colors.pop())
		L.append(g)
		
	plt.xlabel("Expression of gene 1")
	plt.ylabel("Expression of gene 2")
	plt.title("{} ({})".format(GO2T[go], go))
	
	for (scale, loc) in [('linear', 'upper right'), ('log', 'upper left')] :
		plt.xscale(scale)
		plt.yscale(scale)
		
		plt.legend(L, loc=loc)
		
		for ext in PARAM['ext'] :
			plt.savefig(OFILE['2d-plot'].format(go=go_safe, log=scale[0:3], ext=ext))
	
	plt.close()

