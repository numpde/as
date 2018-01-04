
# RA, 2018-01-04

## ================== IMPORTS :

import re
import os
import sys
import math
import pickle
import random
import inspect
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from string import ascii_lowercase
from numpy.matlib import repmat
from scipy import stats
from scipy.constants import golden as phi
from itertools import chain
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress
from networkx.drawing.nx_agraph import graphviz_layout

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'GO=>Info' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'gene-ks'  : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
	'GO=>CI'   : "OUTPUT/0_go2ci/UV/go2ci.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'tree' : "OUTPUT/9_tree/tree_{go}.{ext}",
	'info' : "OUTPUT/9_tree/tree_{go}_info.txt",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	# Top 50 from TXP (2017.01.04)
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
	
	# Last 50 from TXP (2017.01.04)
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
		#"GO:0001525", # angiogenesis
		"GO:0006281", # DNA-repair
		"GO:0006955", # immune response
		#"GO:0007049", # cell cycle
		#"GO:0016477", # cell migration
		#"GO:0004984", # olfactory receptor activity
		"GO:0045596", # negative regulation of cell differentiation
		"GO:0045597", # positive regulation of cell differentiation
	},
	
	# Figure formats
	'ext' : ['png', 'eps', 'pdf'],
	
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.2)),
}

mpl.rcParams['axes.labelsize'] = 'large'

## ====================== AUX :

def abbr(t) :
	D = { 
		"negative regulation" : "-regu",
		"positive regulation" : "+regu",
		"replication" : "repl",
		"regulation" : "regu",
		"involved in" : "in",
		"synthesis" : "syn",
		"double" : "dbl",
		"single" : "sgl",
		"error" : "err",
	}
	
	for (S, s) in sorted(D.items(), key=(lambda x : -len(x[0]))) :
		t = t.replace(S, s)
	
	return t

## ====================== (!) :


## ===================== WORK :

#[ ]#

CI_data = pickle.load(open(IFILE['GO=>CI'], 'rb'))

# GO2E : GO ID --> Clustering index
GO2CI = CI_data['GO2CI']

# GO2E : GO ID --> [ENSG IDs]
GO2E = CI_data['GO2E']

# GO2T : GO ID --> GO category name
GO2T = CI_data['GO2T']



## DEBUG

#G = pickle.load(open(IFILE['GO=>Info'], 'rb'))
#print(G['GO:0031625'])
#exit()
#del G

#[ ]#

class gograph :

	#def edge_type(G, e) :
		#return set(G.get_edge_data(*e).keys())

	#for e in GO.edges() :
		#print(edge_type(GO, e))

	def __init__(self, graph_file) :
		self.G = pickle.load(open(graph_file, 'rb'))
		
		# UV-damage excision repair, DNA incision
		self.G.remove_node('GO:1990731')

		R = {
			# http://amigo.geneontology.org/amigo/term/GO:0008150
			'biological_process',
			
			# http://amigo.geneontology.org/amigo/term/GO:0005575
			'cellular_component',
			
			# http://amigo.geneontology.org/amigo/term/GO:0003674
			'molecular_function',
		}
		#
		R = {
			n : [i for i in self.G.nodes() if (self.G.nodes[i]['name'] == n)]
			for n in R
		}
		#
		assert(all((len(I) == 1) for (n, I) in R.items()))
		#
		self.R = { n : I[0] for (n, I) in R.items() }
		#
		#print("Roots:", self.R)

	def subgraph(self, i, depth=8) :
		I = [i]
		g = nx.MultiDiGraph()
		while (I and (depth > 0) and (g.number_of_nodes() < 1000)) :
			g.add_nodes_from(I)
			E = [(i, j) for i in I for j in sorted(self.G.predecessors(i))]
			E = [(i, j) for (i, j) in E if ("obsolete" not in self.G.nodes[j]['name'])]
			I = [j for (i, j) in E if (j not in g.nodes())]
			g.add_edges_from(E)
			depth -= 1

		g = nx.MultiDiGraph(nx.subgraph(self.G, g.nodes()))
		
		for n in g.nodes() :
			info = g.nodes[n]
			for (k, v) in info.items() :
				if (k == 'name') : continue
				#if (k == 'synonym') : print(k, v)
				g.nodes[n][k] = None
				#print(g.nodes[n])
		
		#exit()
		
		return g


def plot(root, go_filename) :

	print(root, "--", GO2T[root], "({} genes)".format(len(GO2E[root])))
	
	# Ontology subtree
	g = gograph(IFILE['GO=>Info']).subgraph(root)
	
	print(nx.info(g))
	
	#print(root, len(GO2E[root]), GO2E[root])
	
	# For the selection of plt.cm.* colormap, see
	# https://matplotlib.org/examples/color/colormaps_reference.html
	
	node_style = {
		'!' : { 'nodelist' : [], 'node_size' : [], 'node_color' : [] },
		'?' : { 'nodelist' : [], 'node_size' : [], 'node_color' : 'y' },
	}
	
	min_node_size = 5
	
	# Colormap
	cmap = plt.cm.cool
	
	# The ENSG's for each GO term that are also contained in the root GO term
	GO2e = {
		go : GO2E[root] & GO2E[go]
		for go in g.nodes
	}
	
	for go in g.nodes :
		ci = GO2CI[go]
		if ci :
			node_style['!']['nodelist'].append(go)
			node_style['!']['node_size'].append(min_node_size + len(GO2E.get(go, {})))
			node_style['!']['node_color'].append(cmap((ci+1)/2))
		else :
			node_style['?']['nodelist'].append(go)
			node_style['?']['node_size'].append(min_node_size)
	
	pos = graphviz_layout(g, prog='dot', root=root)
	pos = { i : (x, -y) for (i, (x, y)) in pos.items() }
	
	plt.figure(figsize=(10, 15), dpi=150)

	plt.axis('off')
	
	for style in node_style.values() :
		nx.draw_networkx_nodes(g, pos, linewidths=0.1, **style)
	
	nx.draw_networkx_edges(g, pos, width=0.1, arrows=False, edge_color='y')
	
	for (go, (x, y)) in pos.items() :
		al = { 'horizontalalignment' : 'left', 'verticalalignment' : 'bottom' }
		plt.text(x, y, abbr(g.nodes[go]['name'])[0:45] + " ({}, {}/{})".format(go[3:], len(GO2e[go]), len(GO2E[go])), fontsize=3, rotation=60, **al)
	
	plt.title("{} -- {}".format(root, GO2T[root]))
	
	## COLORBAR
	#ax = plt.axes([0.01, 0.01, 0.2, 0.01], facecolor='y')
	#p = list(ax.get_position().bounds)
	#ax.axis('off')
	#ax.imshow(repmat(np.linspace(-1, 1, 256), 2, 1), aspect='auto', cmap=cmap)
	#ax.text(0.01, 0.1, "-1", va='bottom', ha='left',  transform=ax.transAxes, fontsize=7)
	#ax.text(0.99, 0.1, "+1", va='bottom', ha='right', transform=ax.transAxes, fontsize=7)
	#ax.text(0.50, 0.1, "Clustering index", va='bottom', ha='center', transform=ax.transAxes, fontsize=7)
	
	# STATS
	ax = plt.axes([0.7 - 0.01, 0.03, 0.3, 0.1], facecolor='w')
	
	#ax.labelsize = 'small'
	go2size = dict(zip(node_style['!']['nodelist'], node_style['!']['node_size']))
	go2color = dict(zip(node_style['!']['nodelist'], node_style['!']['node_color']))
	for go in g.nodes :
		if (not GO2E[go]) : continue
		plt.semilogx(len(GO2E[go]), GO2CI[go], 'o', color=go2color[go], markersize=math.sqrt(go2size[go]))
	
	plt.ylim(-1, +1) # Possible range of the clustering index
	
	xlim = plt.xlim()
	plt.semilogx(xlim, (GO2CI[root],) * 2, '-k', linewidth=1)
	
	plt.xlabel("Size", fontsize=6)
	plt.ylabel("Clustering index", fontsize=6)
	ax.tick_params(axis='x', labelsize=6)
	ax.tick_params(axis='y', labelsize=6)
	
	# SAVE
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['tree'].format(go=go_filename, ext=ext))
	
	plt.close()



for (root, x) in zip(PARAM['Top50'], ascii_lowercase) :
	go_filename = "Top50_" + x
	plot(root, go_filename)


for (root, x) in zip(PARAM['Last50'], ascii_lowercase) :
	go_filename = "Last50_" + x
	plot(root, go_filename)


for root in sorted(PARAM['GO filter']) :
	go_filename = str(root).replace(':', '-')
	plot(root, go_filename)

