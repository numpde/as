
# RA, 2018-01-04

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

from collections     import defaultdict
from string          import ascii_lowercase
from numpy.matlib    import repmat
from scipy           import stats
from scipy.constants import golden as phi
from itertools       import chain
from multiprocessing import cpu_count
from joblib          import Parallel, delayed
from progressbar     import ProgressBar as Progress

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
	#'tree' : "OUTPUT/9_tree/tree_{go}.{ext}",
	#'info' : "OUTPUT/9_tree/tree_{go}_info.txt",
	
	# Overview figure
	'cisz' : "OUTPUT/B_go-by-wq/ci-vs-sz_{p}.{ext}",
	
	'list lo' : "OUTPUT/B_go-by-wq/list_lo.txt",
	'list hi' : "OUTPUT/B_go-by-wq/list_hi.txt",
	
	# All GO terms 
	'list all' : "OUTPUT/B_go-by-wq/list_all.txt",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	# Top 50 from TXP (2017-01-04)
	'Top50' : [
		"GO:0000922", # spindle pole
		"GO:0007062", # sister chromatid cohesion
		"GO:0006260", # DNA replication
		"GO:0012505", # endomembrane system
		"GO:0005814", # centriole
		"GO:0005085", # guanyl-nucleotide exchange factor activity
		"GO:0018105", # peptidyl-serine phosphorylation
		"GO:0015630", # microtubule cytoskeleton
		"GO:0017137", # Rab GTPase binding
		"GO:0005815", # microtubule organizing center
	],
	
	# Last 50 from TXP (2017-01-04)
	'Last50' : [
		"GO:0007275", # multicellular organism development
		"GO:0016567", # protein ubiquitination
		"GO:0019901", # protein kinase binding
		"GO:0016032", # viral process
		"GO:0045296", # cadherin binding
		"GO:0006357", # regulation of transcription from RNA polymerase II promoter
		"GO:0005622", # intracellular
		"GO:0030154", # cell differentiation
		"GO:0031625", # ubiquitin protein ligase binding
		"GO:0046982", # protein heterodimerization activity
		"GO:0043234", # protein complex
		"GO:0007155", # cell adhesion
		"GO:0008284", # positive regulation of cell proliferation
		"GO:0045892", # negative regulation of transcription, DNA-templated
		"GO:0043066", # negative regulation of apoptotic process
		"GO:0019899", # enzyme binding
		"GO:0008283", # cell proliferation
		"GO:0005743", # mitochondrial inner membrane
		"GO:0043312", # neutrophil degranulation
		"GO:0042493", # response to drug
	],
	
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
	
	# Quantiles for extrema
	'q lo' : 0.05,
	'q hi' : 0.95,
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

pass

## ===================== DATA :

#[ ]#

# Clustering indices data bundle
CI_data = pickle.load(open(IFILE['GO=>CI'], 'rb'))

# GO2E : GO ID --> Clustering index
GO2CI = CI_data['GO2CI']

# GO2E : GO ID --> [ENSG IDs]
GO2E = CI_data['GO2E']

# GO2T : GO ID --> GO category name
GO2T = CI_data['GO2T']

# N2CI : size of GO term --> [clustering indices]
N2CI = CI_data['N2CI']

# GO2WQ : GO ID --> windowed quantile
GO2WQ = CI_data['GO2WQ']

## The Gene Ontology graph
#GO_graph = pickle.load(open(IFILE['GO graph'], 'rb'))

## Are those GO IDs in the GO graph?
#go_not_in_graph = set(GO2E.keys()) - set(GO_graph.nodes())
#print("Note: {} GO IDs are not in the graph".format(len(go_not_in_graph)))

# GO categories that have windowed quantiles
GO_all = [go for (go, q) in GO2WQ.items() if q]
# Order them by size
GO_all = sorted(GO_all, key=(lambda go : len(GO2E[go])))

# Extrema and bulk
GO_lo = [ go for go in GO_all if (GO2WQ[go] <= PARAM['q lo']) ]
GO_hi = [ go for go in GO_all if (GO2WQ[go] >= PARAM['q hi']) ]
GO_in = list(set(GO_all) - (set(GO_lo) & set(GO_hi)))

## ===================== WORK :

#[ ]#


def plot_overview(t=None) :
	
	# t = threshold for percentile
	if (t is None) :
		for t in [0, 1, 2, 5, 10, 90, 95, 98, 99] :
			plot_overview(t)
		return
	
	plt.close('all')

	plt.figure(figsize=(math.ceil(6*phi), 6), dpi=150)
	
	cmap = plt.cm.cool

	P = [(len(GO2E[go]), GO2CI[go]) for go in GO_all if (GO2WQ[go] > t/100)]
	plt.plot(*zip(*P), '.', color='r', markersize=3, zorder=0)
	
	if t :
		P = [(len(GO2E[go]), GO2CI[go]) for go in GO_all if (GO2WQ[go] <= t/100)]
		plt.plot(*zip(*P), 'x', color='b', markersize=3, zorder=1)
		
		plt.legend(["Above {}%".format(t), "Below {}%".format(t)])
	
	#plt.scatter(*zip(*[(len(GO2WQ[go]), GO2CI[go]) for go in GO_all]), c=[GO2WQ[go] for go in GO_all], cmap=cmap)
	
	plt.xscale('log')
	
	plt.ylim((-1, 1))
	
	plt.xlabel("Size of GO category")
	plt.ylabel("Clustering index")
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['cisz'].format(p=t, ext=ext))
	
	#plt.show()
	plt.close('all')


def dump_extrema() :
	
	def write(GO, f, header=True) :
		if header :
			print("GO ID", "GO name", "Size", "Clustering index", "CI quantile", sep='\t', file=f) 
		for go in GO :
			print(go, GO2T[go], len(GO2E[go]), GO2CI[go], GO2WQ[go], sep='\t', file=f)
	
	write(GO_lo, open(OFILE['list lo'], 'w'))
	write(GO_hi, open(OFILE['list hi'], 'w'))

	write(sorted(GO_all, key=(lambda go : GO2WQ[go])), open(OFILE['list all'], 'w'))

###

def main() :
	plot_overview()
	dump_extrema()


main()
