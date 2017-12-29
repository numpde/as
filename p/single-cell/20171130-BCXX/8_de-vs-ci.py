
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
	'2d-plot' : "OUTPUT/8_de-vs-ci/de-vs-ex_{go}.{ext}",
	'2d-info'  : "OUTPUT/8_de-vs-ci/de-vs-ex_{go}_info.txt",
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
}

mpl.rcParams['axes.labelsize'] = 'large'

## ====================== (!) :


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
S = [ np.array(G2S[g]) for g in G ]
del G2S # Use G and S


#[ COMPUTE EXPRESSION MAGNITUDE ]#

E2EX = { e : np.mean(E2X[e]) for e in BC_E }


#[ LOAD GO TERMS & ANNOTATION ]#

# GO2E : GO ID --> [ENSG IDs]
# Read it from file
GO2E = {
	go_E[0] : go_E[1:]
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


#[ FOR EACH GO TERM OF INTEREST... ]#

for go in PARAM['GO filter'] :
	
	# GO ID for file names
	go_safe = str(go).replace(':', '-')


	#[ SELECT GENES FOR THIS GO TERM ]#

	# Default choice: all genes
	E = set(BC_E)
	
	# Filter by GO
	if (go is not None) : E &= set(GO2E[go])
	
	# Order by differential expression
	E = list(sorted(E, key=(lambda e : -E2DE[e])))
	
	# Pick the top 2 for a 2d plot
	E = E[0:2]
	
	# Get the indices of those genes
	E = [BC_E.index(e) for e in E]
	
	#[ PLOT 2-GENE EXPRESSION ]#
	
	plt.figure(figsize=(math.ceil(6*phi), 6), dpi=150)
	
	colors = plt.get_cmap('hsv')(np.linspace(0.1, 1.0, len(G))).tolist()
	
	L = []
	for (g, s) in zip(G, S) :
		Y = np.take(X, E, axis=axis_gene)
		Y = np.take(Y, s, axis=axis_smpl)
		a = np.take(Y, 0, axis=axis_gene)
		b = np.take(Y, 1, axis=axis_gene)
		plt.loglog(a, b, '.', markersize=8, color=colors.pop())
		L.append(g)
	
	plt.legend(L, loc='upper left')
	plt.xlabel("Expression gene 1")
	plt.ylabel("Expression gene 2")
	plt.title("{} ({})".format(GO2T[go], go))
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['2d-plot'].format(go=go_safe, ext=ext))
	
	plt.close()
