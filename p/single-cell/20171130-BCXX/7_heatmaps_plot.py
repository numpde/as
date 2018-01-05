
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

from scipy           import stats
from itertools       import chain
from scipy.constants import golden      as phi
from progressbar     import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'GO=>Info' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'gene-ks'  : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
}

## =================== OUTPUT :

OFILE = {
	'de-vs-ex' : "OUTPUT/7_heatmaps/de-vs-ex_{go}.{ext}",
	'dx-info'  : "OUTPUT/7_heatmaps/de-vs-ex_{go}_info.txt",
	
	'heatmap'  : "OUTPUT/7_heatmaps/heatmap_{go}.{ext}",
	'hm-info'  : "OUTPUT/7_heatmaps/heatmap_{go}_info.txt",
	
	'ex-sect'  : "OUTPUT/7_heatmaps/ex-sect_{n}.{ext}",
	'xs-info'  : "OUTPUT/7_heatmaps/ex-sect_{n}_info.txt",
}

## ==================== PARAM :

PARAM = {
	# GO terms of interest
	'GO filter' : {
		None, 
		"GO:0032201", # telomere maintenance via semi-conservative replication
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
		"GO:0043065", # positive regulation of apoptotic process
		"GO:0043066", # negative regulation of apoptotic process
		"GO:0007569", # cell aging
		"GO:0006915", # apoptotic process
		"GO:0012501", # programmed cell death
		"GO:0070265", # necrotic cell death
		"GO:0031724", # CXCR5 chemokine receptor binding
		"GO:0002039", # p53 binding
		"GO:0030330", # DNA damage response, signal transduction by p53 class mediator
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


#[ FOR EACH SECTION OF EXPRESSION MEAN... #]


for quantile in ['-', '+'] :
	
	ab = np.logspace(
		math.log10(min(E2EX.values())) - 0.01, 
		math.log10(max(E2EX.values())) + 0.01, 
		20
	)
	
	for (n, (a, b)) in enumerate(zip(ab[:-1], ab[1:])) :
		
		plt.figure(figsize=(math.ceil(8*phi), 8), dpi=150)
		
		# Background
		if True :
			de = [ E2DE[e] for e in BC_E ]
			ex = [ E2EX[e] for e in BC_E ]
			plt.loglog(ex, de, 'r.', markersize=0.5, zorder=-20)
		
		# Selection of genes in the expression window
		E = { e for e in BC_E if (a <= E2EX[e] < b) }
		
		if (not E) : continue

		q = 5
		pa = np.percentile([E2DE[e] for e in E], q)
		pb = np.percentile([E2DE[e] for e in E], 100-q)
		
		# Box
		if True : 
			xlim = plt.gca().get_xlim()
			ylim = plt.gca().get_ylim()
			plt.loglog([a, b, b, a, a], [pa, pa, pb, pb, pa], '-k', zorder=-30)
			plt.loglog([a, a], ylim, '-k', zorder=-30)
			plt.loglog([b, b], ylim, '-k', zorder=-30)
			plt.gca().set_xlim(xlim)
			plt.gca().set_ylim(ylim)

		assert(quantile in ['-', '+'])
		if (quantile == '-') : E = { e for e in E if (E2DE[e] <= pa) }
		if (quantile == '+') : E = { e for e in E if (E2DE[e] >= pb) }
		
		GO2F = { go : (set(F) & E)  for (go, F) in GO2E.items() }
		GO2F = { go : F for (go, F) in GO2F.items() if F }
		# Sort by most specific GO terms (i.e. slimmest)
		GO2F = list(sorted(
			GO2F.items(), 
			key=(lambda x : (len(GO2E[x[0]]), x[0]))
		))
		
		H = [] # Handles
		L = [] # Legend
		
		# Genes with known GO ID
		for (go, F) in GO2F :
			de = [ E2DE[e] for e in F ]
			ex = [ E2EX[e] for e in F ]
			h2 = plt.loglog(ex, de, 'bo', markersize=4)[0]
			H.append(h2)
			L.append(GO2T.get(go, go)[0:50] + " ({}/{})".format(len(F), len(GO2E[go])))
		
		# If too many legend items...
		max_leg_len = 30
		if (len(L) > max_leg_len) :
			H = H[0:(max_leg_len+1)]
			L = L[0:max_leg_len] + ["(And {} more GO terms...)".format(len(L) - max_leg_len)]
		
		# Genes without GO ID
		U = set(E) - set(chain.from_iterable(dict(GO2F).values()))
		if U :
			de = [ E2DE[e] for e in U ]
			ex = [ E2EX[e] for e in U ]
			h3 = plt.loglog(ex, de, 'rx', markersize=4, zorder=-10)[0]
			H.append(h3)
			L.append("No GO ID ({} genes)".format(len(U)))
		
		if H :
			plt.legend(H, L, loc='upper left', prop={'size': 8})
		
		plt.xlabel('Overall expression mean')
		plt.ylabel('Inter-tumor differential expression (KS {})'.format(KS_meta))
		
		for ext in PARAM['ext'] :
			plt.savefig(OFILE['ex-sect'].format(n=(str(n)+quantile), ext=ext))
		
		plt.close()


#[ FOR EACH GO TERM OF INTEREST... ]#

for go in PARAM['GO filter'] :
	
	# GO ID for file names
	go_safe = str(go).replace(':', '-')


	#[ SELECT GENES FOR THIS GO TERM ]#

	# Default choice: all genes
	E = set(BC_E)
	
	# Filter by GO
	if (go is not None) : E &= set(GO2E[go])
	
	
	#[ PLOT EXPRESSION VS DE ]#
	
	plt.figure(figsize=(math.ceil(8*phi), 8), dpi=150)

	# Background
	de = [ E2DE[e] for e in BC_E ]
	ex = [ E2EX[e] for e in BC_E ]
	plt.loglog(ex, de, 'r.', markersize=0.5)

	if go : 
		de = [ E2DE[e] for e in E ]
		ex = [ E2EX[e] for e in E ]
		h2 = plt.loglog(ex, de, 'bo', markersize=4)[0]
		plt.legend([h2], ["{} ({})".format(GO2T[go], len(E))], prop={'size': 12})
	
	
	plt.xlabel('Overall expression mean')
	plt.ylabel('Inter-tumor differential expression (KS {})'.format(KS_meta))
	
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

	plt.close()
	

	#[ PLOT HEATMAP ]#
	

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
	
	with open(OFILE['hm-info'].format(go=go_safe), 'w') as f :
		print("{} ENSG IDs in selected BC data".format(len(E)), file=f)
		print("GO filter: {}".format(go), file=f)
		if (go is None) : continue
		print(GO2T[go], file=f)
		print("{} ENSG IDs".format(len(GO2E[go])), file=f)

