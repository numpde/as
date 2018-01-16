
# RA, 2018-01-15

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

# http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
from sklearn.manifold import TSNE

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
from scipy.cluster.hierarchy import linkage
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.leaves_list.html
from scipy.cluster.hierarchy import leaves_list

from sklearn.metrics.pairwise import euclidean_distances as euc_dist

# 
from networkx.drawing.nx_agraph import graphviz_layout

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO graph' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'GO=>CI'   : "OUTPUT/0_go2ci/UV/go2ci.pkl",
	
	'tsne runs'  : "OUTPUT/C_goordinates/tsne_runs.pkl",
	'classified' : "OUTPUT/D_classifier-nn/classified.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	
	'tsne'      : "OUTPUT/C_goordinates-D/tsne_dim={dim}_run={run}_sub={sub}.{ext}",
	'tsne info' : "OUTPUT/C_goordinates-D/tsne_dim={dim}.txt",
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

# Abbreviate a GO term t, truncating it to max_len characters
def abbr(t, max_len) :
	D = { 
		"replication" : "repl",
		"involved in" : "in",
		"regulation" : "regu",
		"synthesis" : "syn",
		"negative" : "neg",
		"positive" : "pos",
		"double" : "dbl",
		"single" : "sgl",
		"error" : "err",
	}
	
	for pair in D.items() : t = t.replace(*pair)

	if (len(t) > max_len) : t = t[0:(max_len-3)] + "..."
	
	return t

# This script
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

## ====================== (!) :

def goordinate_trafo(GO_I, n_genes) :
	assert(type(GO_I) is list)
	
	T = lil_matrix((len(GO_I), n_genes))
	
	for (n, (go, I)) in enumerate(GO_I) : 
		T[n, I] = 1
	
	return T


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
X = BC_data['X']

# Rearrange data axes
(axis_gene, axis_smpl) = (0, 1)
X = np.moveaxis(X, BC_data['axis_gene'], axis_gene)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# ENSG IDs
BC_E = BC_data['gene_id']
assert(len(BC_E) == X.shape[axis_gene]), "Inconsistent gene info"

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


## The Gene Ontology graph
#GO_graph = pickle.load(open(IFILE['GO graph'], 'rb'))

## Are those GO IDs in the GO graph?
#go_not_in_graph = set(GO2E.keys()) - set(GO_graph.nodes())
#print("Note: {} GO IDs are not in the graph".format(len(go_not_in_graph)))


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

#[ ]#


def data_in_go_space(N=None, norm_features=True, norm_samples=True) :
	
	# Move to the GO feature space
	assert(axis_gene == 0)
	Y = goordinate_trafo(sorted(GO2I.items()), n_genes) * X
	
	# Nontrivial features
	Y = [(go, x) for (go, x) in zip(sorted(GO2I.keys()), Y.tolist()) if GO2WQ.get(go, None) and np.sum(x)]
	# Sort by "clustering index windowed quantile"
	Y = sorted(Y, key=(lambda go_x : GO2WQ[go_x[0]]))
	
	# Split GO ID / Data 
	(GO, Y) = zip(*Y)
	
	# Filter by size
	(GO, Y) = zip(*[(go, y) for (go, y) in zip(GO, Y) if (len(GO2E[go]) <= 9)])
	
	# Take the first N 
	if N :
		(GO, Y) = zip(*[(go, y) for (go, y) in list(zip(GO, Y))[0:N]])
	
	# Collect the features into a numpy matrix row-wise
	Y = np.vstack(Y)
	
	# Normalize data feature-wise
	if norm_features : 
		Y = np.vstack((y / (np.sum(y) or 1)) for y in Y)
	
	## Normalize data sample-wise
	#if norm_samples :
		#Y = np.vstack((s / (np.sum(s) or 1)) for s in Y.transpose()).transpose()
	
	return (GO, Y)


def plot_CI_vs_N() :
	pass



def compute_all_tSNE_in_go_space() :
	
	runs_filename = IFILE['tsne runs']
	
	assert(os.path.isfile(runs_filename)), ("File {} with t-SNE not found.".format(runs_filename))

	return pickle.load(open(runs_filename, 'rb'))['runs']


def plot_tSNE_in_go_space() :
	
	# Predicted class
	Yp = pickle.load(open(IFILE['classified'], 'rb'))['Yp']
	assert(Yp.shape[0] == n_samples)
	
	for run_info in compute_all_tSNE_in_go_space() :
		
		Z   = run_info['Z']
		N   = run_info['N']
		run = run_info['run']
		#     run_info['GO']
	
		# Log the selected GO terms to file (including redundancies)
		with open(OFILE['tsne info'].format(dim=N), 'w') as f :
			for go in run_info['GO'] :
				print(go, len(GO2E[go]), GO2T[go], " / ".join([(i + " -- " + GO2T[i]) for i in GO2A[go][1:]]), sep='\t', file=f)
		
		cm_a = plt.cm.winter
		cm_b = plt.cm.cool
		cm_c = plt.cm.autumn_r
		
		# Legend and colors for the tumors
		(L, c) = zip(*[
			( "BC01 (ER+)",         cm_a(0.0) ),
			( "BC02 (ER+)",         cm_a(0.4) ),
			
			( "BC03 (ER+, HER2+)",  cm_a(0.7) ),
			( "BC03LN",             cm_a(1.0) ),
			
			( "BC04 (HER2+)",       cm_b(0.3) ),
			( "BC05 (HER2+)",       cm_b(0.5) ),
			( "BC06 (HER2+)",       cm_b(0.9) ),
			
			( "BC07 (TNBC)",        cm_c(0.0) ),
			( "BC07LN",             cm_c(0.1) ),
			( "BC08 (TNBC)",        cm_c(0.3) ),
			( "BC09 (TNBC)",        cm_c(0.5) ),
			( "BC09_Re",            cm_c(0.6) ),
			( "BC10 (TNBC)",        cm_c(0.8) ),
			( "BC11 (TNBC)",        cm_c(1.0) ),
		])
		#
		# Check that the legend corresponds to the cell groups
		assert(all((l.startswith(k) for (l, k) in zip(L, sorted(G2S.keys())))))
		#
		# The number of samples will be filled in later
		L = [("{} x " + l) for l in L]
		
		plt.close('all')
		
		# Clean the axes
		plt.xticks([], [])#; plt.xlabel("t-SNE 1")
		plt.yticks([], [])#; plt.ylabel("t-SNE 2")
		
		# Handles and texts for the legend
		HL = []
		
		for background in [True, False] :
			# Iterate over the cell groups
			for (n, (g, s)) in enumerate(sorted(G2S.items())) :
				
				# NN classification of this sample into this group
				sz = 1 + (30 * Yp[s, n])
			
				if background :
					# Scatter plot with no color
					plt.scatter(*Z[:, s], s=sz, facecolors='None', edgecolors='k', lw=0.1)
				
				else :
					# For one particular configuration:
					# Save intermediate plots showing individual groups better
					if ( (N, run) == (20, 1) ) :
						for ext in PARAM['ext'] : 
							plt.savefig(OFILE['tsne'].format(dim=N, run=run, sub=n, ext=ext))
					
					# Keep only the "healthy" cells, assuming most cells are healthy
					s = [c for c in s if (np.mean(X[:, c]) >= np.median(X[:, c])/2)]
					
					# Scatter plot of the cell group
					h = plt.scatter(*Z[:, s], alpha=0.8, c=c[n], s=sz, edgecolors='k', lw=0.2)
					
					# Fill in the number of samples in the legend
					HL.append(( h, L[n].format(len(s)) ))
					
					plt.legend(*zip(*HL), prop={'size': 5}, loc='upper left')
		
		for ext in PARAM['ext'] : 
			plt.savefig(OFILE['tsne'].format(dim=N, run=run, sub="all", ext=ext))

###

if (__name__ == "__main__") :
	plot_tSNE_in_go_space()
