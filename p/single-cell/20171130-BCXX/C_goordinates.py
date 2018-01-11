
# RA, 2018-01-11

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

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
from scipy.cluster.hierarchy import linkage
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.leaves_list.html
from scipy.cluster.hierarchy import leaves_list

# 
from networkx.drawing.nx_agraph import graphviz_layout

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO graph' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'GO=>CI'   : "OUTPUT/0_go2ci/UV/go2ci.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'tsne' : "OUTPUT/C_goordinates/tsne_dim={dim}_run={run}.{ext}",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	
	# Further GO terms of interest
	'GO filter' : {
		"GO:0001525", # angiogenesis
		"GO:0006281", # DNA-repair
		"GO:0006955", # immune response
		"GO:0007049", # cell cycle
		"GO:0016477", # cell migration
		"GO:0004984", # olfactory receptor activity
		"GO:0045596", # negative regulation of cell differentiation
		"GO:0045597", # positive regulation of cell differentiation
		"GO:0000723", # telomere maintenance
		
		"GO:0007173", # epidermal growth factor receptor signaling pathway
		"GO:0035004", # phosphatidylinositol 3-kinase activity
		"GO:0051019", # mitogen-activated protein kinase binding
	},
	
	# Figure formats
	'ext' : ['png', 'pdf'],
	
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.2)),
}

mpl.rcParams['axes.labelsize'] = 'large'

## ====================== AUX :

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

## ====================== (!) :

def goordinate_trafo(GO_I, n_genes) :
	assert(type(GO_I) is list)
	
	T = lil_matrix((len(GO_I), n_genes))
	
	for (n, (go, I)) in enumerate(GO_I) : T[n, I] = 1
	
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
def cos_dist(X, axis) :
	# Covariance & norm products
	C = np.tensordot(X, X, axes=([axis], [axis]))
	V = np.sqrt(np.outer(np.diag(C), np.diag(C)))
	V[V == 0] = 1
	D = 1 - (C / V)
	return D

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

H2GO = defaultdict(set)
for (go, E) in GO2E.items() : H2GO[hash('+'.join(sorted(E)))].add(go)
H2GO = dict(H2GO)
#
# Check for no hash collisions
assert(all((1 == len(set('+'.join(sorted(GO2E[go])) for go in GO))) for GO in H2GO.values()))
#
# Non-redundant GO categories and their aliases
GO2A = { min(GO) : sorted(GO) for GO in H2GO.values() }
#print("{} of {} GO categories are non-redundant".format(len(GO2A), len(GO2E)))
#
del H2GO

def restrict(GO2X, GO) :
	return { go : x for (go, x) in GO2X.items() if (go in GO) }

GO2E  = restrict(GO2E, GO2A.keys())
GO2I  = restrict(GO2I, GO2A.keys())
GO2CI = restrict(GO2CI, GO2A.keys())
GO2WQ = restrict(GO2WQ, GO2A.keys())

print(len(GO2E))

## ===================== WORK :

#[ ]#


def CI(X) :
	return np.mean([
		np.sign(x) 
		for x in chain.from_iterable(silhouette(cos_dist(X, axis_gene), S).values())
	])

def go_space() :
	
	# Move to the GO feature space
	TX = goordinate_trafo(sorted(GO2I.items()), n_genes) * X
	
	# 
	Y = np.asarray(TX)
	
	# Cut off values <= 1
	#Y = np.zeros(TX.shape)
	#Y[TX > 1] = TX[TX > 1]
	
	# Log-transform
	#Y = np.log(Y)
	
	# Nontrivial features
	Y = [(go, x) for (go, x) in zip(sorted(GO2I.keys()), Y.tolist()) if GO2WQ.get(go, None) and np.sum(x)]
	# Sort by "clustering index windowed quantile"
	Y = sorted(Y, key=(lambda go_x : GO2WQ[go_x[0]]))
	# Split GO ID / Data 
	(GO, Y) = zip(*Y)
	
	return (GO, Y)

def plot_ci2() :
	(GO, Y) = go_space()
	
	
	# Normalize data feature-wise
	Y = [x/np.sum(x) for x in Y]
	
	## Z-normalize data feature-wise
	#Y = [stats.mstats.zscore(x) for x in Y]
	
	plt.ion()
	plt.show()
	n = 1
	while (n < len(GO)) :
		ci = CI(np.vstack(Y[0:math.ceil(n)]))
		plt.semilogx(n, ci, 'o', c=('b' if (ci < 0) else 'r'))
		plt.pause(0.01)
		n *= math.sqrt(2)

def plot_go_space(N) :
	(GO, Y) = go_space()
	
	# Filter GO categories by size
	(GO, Y) = zip(*[(go, y) for (go, y) in zip(GO, Y) if (1 == len(GO2E[go]))])
	
	(GO, Y) = zip(*[(go, y) for (go, y) in list(zip(GO, Y))[0:N]])
	
	#for go in GO :
		#print(go, "({})".format(len(GO2E[go])), GO2T[go])
		#for a in GO2A[go] :
			#if (a != go) :
				#print("aka: ", a, GO2T[a])
	
	# Normalize data feature-wise
	Y = [y/np.sum(y) for y in Y]
	
	Y = np.vstack(Y)
	#Y[Y != 0] = 1
	
	# Normalize data sample-wise
	Y = np.vstack(s / (np.sum(s) or 1) for s in Y.transpose()).transpose()
	
	## Iterate over samples/cells
	#for (n, f) in enumerate(Y.transpose().tolist()) :
		## Pick the top N features
		#N = 3
		#Y[(Y[:, n] > sorted(f, reverse=True)[N]), n] = 1
	##
	#Y[Y < 1] = 0
	
	#for (n, f) in enumerate(Y.transpose().tolist()) :
		#print("Cell ", n)
		#for i in np.nonzero(Y[:, n])[0] : 
			#print(GO2T[GO[i]])
	
	c = [
		plt.cm.winter(0.0),  # BC01:    ER+
		plt.cm.winter(0.4),  # BC02:    ER+
		
		plt.cm.winter(0.8),  # BC03:    ER+ and HER2+
		plt.cm.winter(1.0),  # BC03LN
		
		plt.cm.cool  (0.3),  # BC04:    HER2+
		plt.cm.cool  (0.5),  # BC05:    HER2+
		plt.cm.cool  (0.9),  # BC06:    HER2+
		
		plt.cm.autumn(0.0),  # BC07:    TNBC
		plt.cm.autumn(0.1),  # BC07LN
		plt.cm.autumn(0.3),  # BC08:    TNBC
		plt.cm.autumn(0.5),  # BC09:    TNBC
		plt.cm.autumn(0.6),  # BC09_Re
		plt.cm.autumn(0.8),  # BC10:    TNBC
		plt.cm.autumn(1.0),  # BC11:    TNBC
	]
	
	
	#print(X.shape)
	#for (n, (g, s)) in enumerate(sorted(G2S.items())) :
		#print(n, [np.mean(X[:, c]) for c in s])
	
	
	for run in range(5) : 
		plt.close('all')
		
		from sklearn.manifold import TSNE
		Z = TSNE(n_components=2).fit_transform(Y.transpose()).transpose()
		#
		for (n, (g, s)) in enumerate(sorted(G2S.items())) :
			s = [c for c in s if (np.mean(X[:, c]) >= 10)]
			plt.scatter(*Z[:, s], alpha=0.8, c=c[n], s=7) # c=plt.cm.tab20(n/19)
		#
		plt.legend(list(sorted(G2S.keys())), prop={'size': 6}, loc='upper left')
		plt.axis('off')
		
		for ext in PARAM['ext'] : 
			plt.savefig(OFILE['tsne'].format(dim=N, run=run, ext=ext))
		
		#print(n)
		#for t in [GO2T[go] for (_, go) in sorted(zip(f, GO), reverse=True)[0:5]] :
			#print(t)
	
	## Reorder the features by similarity
	#I = list(leaves_list(linkage(Y, method='complete')))
	#GO = [GO[i] for i in I]
	#Y = np.vstack([Y[i] for i in I])
	
	## Lexicographic sort by feature
	#Y = Y[np.lexsort(np.rot90(Y))]
	
	
	## Reorder samples
	#Y = Y[:, list(leaves_list(linkage(Y.transpose(), method='complete')))]
	
	#plt.imshow(Y, aspect='auto')
	#plt.show()

###

def main() :
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(plot_go_space)(N)
		for N in [10, 15, 20, 30, 40, 60, 100, 200]
	)
	
	
	#plot_ci2()
	#input()
	


main()
