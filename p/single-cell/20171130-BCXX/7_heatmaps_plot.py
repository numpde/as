
# RA, 2017-12-20

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy.constants import golden as phi
from itertools import chain
from progressbar import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'GO=>Info' : "ORIGINALS/go/go-summary.csv",
	'gene-ks'  : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
}

## =================== OUTPUT :

OFILE = {
	'de-vs-ex' : "OUTPUT/7_heatmaps/de-vs-ex_{go}.{ext}",
	'heatmap'  : "OUTPUT/7_heatmaps/heatmap_{go}.{ext}",
	'go-info'  : "OUTPUT/7_heatmaps/info_{go}.txt",
}

## ==================== PARAM :

PARAM = {
	# GO terms of interest
	'GO filter' : {
		None, 
		"GO:0001525", # angiogenesis
		"GO:0006281", # DNA-repair
		"GO:0006955", # immune response
		"GO:0007049", # cell cycle
		"GO:0016477", # cell migration
		"GO:0004984", # olfactory receptor activity
		"GO:0004930", # G-protein coupled receptor activity
		"GO:0070531", # BRCA1-A complex
		"GO:0070532", # BRCA1-B complex
		"GO:0070533", # BRCA1-C complex
		"GO:0006915", # apoptotic process
		"GO:0043065", # positive regulation of apoptotic process
		"GO:0043066", # negative regulation of apoptotic process
	},
	
	# Figure formats
	'ext' : ['png'],
}

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

# Batch to (sample index, header)
B2SH = BC_data['B2SH']

# Clusters/groups by batch
G = sorted(list(B2SH.keys()))
S = [ np.array([s for (s, h) in B2SH[g]]) for g in G ]


#[ LOAD KS DIFFERENTIAL EXPRESSION DATA ]#

KS_data = pickle.load(open(IFILE['gene-ks'], 'rb'))
#print(KS.keys())

KS_OP = 1
KS_meta = KS_data['KS_meta'][KS_OP]
assert(KS_meta in ['max', 'mean', 'median', 'min', 'std'])

E2DE = { e : ks[KS_OP] for (e, ks) in KS_data['E2KS'].items() }


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
GO2T = dict(
	tuple(L.rstrip().split('\t')[:2])
	for L in open(IFILE['GO=>Info'], 'r').readlines()[1:]
)

for go in list(PARAM['GO filter']) :
	if (go is not None) and (go not in GO2E) :
		print("NOTE: {} has no gene list".format(go))
		PARAM['GO filter'].remove(go)

#[ FOR EACH GO TERM OF INTEREST ]#

for go in PARAM['GO filter'] :
	
	# GO ID for file names
	go_safe = str(go).replace(':', '-')


	#[ SELECT GENES OF INTEREST ]#

	# Default choice: all genes
	E = set(BC_E)
	
	# Filter by GO
	if (go is not None) : E &= set(GO2E[go])
	
	
	#[ PLOT EXPRESSION VS DE ]#
	
	plt.clf()

	# E2X : BC ENSG --> BC X data
	E2X = dict(zip(BC_E, np.moveaxis(X, axis_gene, 0)))

	# Background
	de = [ E2DE[e]  for e in BC_E ]
	ex = [ np.mean(E2X[e]) for e in BC_E ]
	plt.loglog(ex, de, 'r.', markersize=0.4)

	if go : 
		de = [ E2DE[e]  for e in E ]
		ex = [ np.mean(E2X[e]) for e in E ]
		h2 = plt.loglog(ex, de, 'bo', markersize=2)[0]
		plt.legend([h2], ["{} ({})".format(GO2T[go], len(E))])
	
	
	plt.xlabel('Expression mean')
	plt.ylabel('Differential expression (KS {})'.format(KS_meta))
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['de-vs-ex'].format(go=go_safe, ext=ext))
	
	
	#[ SELECT GENES FOR HEATMAP ]#
	
	# Order by differential expression
	E = list(sorted(E, key=(lambda e : -E2DE[e])))
	
	print("GO filter: {}, # genes: {}".format(go, len(E)))
	
	Y = np.take(X, [BC_E.index(e) for e in E], axis=axis_gene)
	
	
	#[ IN-CLUSTER REORDER FOR IMPROVED VISUALS ]#

	assert(axis_smpl == 0), "General case not implemented"
	for s in S :
		Z = Y[s, :]
		Z = Z[sorted(range(Z.shape[0]), key=(lambda i : -np.sum(Z[i, :]))), :]
		Y[s, :] = Z


	#[ PREPARE HEATMAP MATRIX ]#

	Z = np.zeros(Y.shape)
	#
	Z[Y != 0] = np.log(Y[Y != 0])
	Z[Z <= 0] = 0
	Z /= np.max(Z)
	#
	Z[Y == 0] = -1
	#
	Y = Z


	#[ PLOT HEATMAP ]#
	
	plt.clf()

	# 
	R = [len(s) for s in S]
	spacing = 0
	R = [0] + list(chain.from_iterable([[spacing, b] for b in R]))[1:]
	R = np.cumsum(R)
	R = R / max(R)
	
	fig = plt.figure(figsize=(math.ceil(4*phi), 4), dpi=300)

	AX = [
		fig.add_axes([a, 0, b-a, 1])
		for (a, b) in zip(R[0::2], R[1::2])
	]

	for (n, ax) in enumerate(AX) :
		y = np.take(Y, S[n], axis=axis_smpl)
		ax.matshow(
			np.moveaxis(y, axis_gene, 0), 
			aspect='auto', origin='upper', cmap="YlGnBu", 
			vmin=-0.1, vmax=1
		)
		
		ax.set_xticks([])
		ax.set_yticks([])
		
		# Group label
		ax.text(
			0.95, 0.005, G[n], 
			size=6,
			transform=ax.transAxes, rotation=90, 
			horizontalalignment='right', verticalalignment='bottom',
		)
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['heatmap'].format(go=go_safe, ext=ext))
	
	plt.close()
	
	
	#[ WRITE INFO ABOUT THE GO TERM ]#
	
	with open(OFILE['go-info'].format(go=go_safe), 'w') as f :
		print("{} ENSG IDs in selected BC data".format(len(E)), file=f)
		print("GO filter: {}".format(go), file=f)
		if (go is None) : continue
		print(GO2T[go], file=f)
		print("{} ENSG IDs".format(len(GO2E[go])), file=f)

