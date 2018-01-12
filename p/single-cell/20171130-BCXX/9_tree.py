
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
	'tree' : "OUTPUT/9_tree/tree_{go}.{ext}",
	'info' : "OUTPUT/9_tree/tree_{go}_info.txt",
	'cisz' : "OUTPUT/9_tree/ci-vs-sz.{ext}",
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
}

mpl.rcParams['axes.labelsize'] = 'large'

## ====================== AUX :

def abbr(t, n) :
	D = { 
		"negative regulation" : "neg regu",
		"positive regulation" : "pos regu",
		"replication" : "repl",
		"involved in" : "in",
		"regulation" : "regu",
		"synthesis" : "syn",
		"double" : "dbl",
		"single" : "sgl",
		"error" : "err",
	}
	
	for (S, s) in sorted(D.items(), key=(lambda x : -len(x[0]))) :
		t = t.replace(S, s)
	
	if (len(t) > (n + 3)) : 
		t = t[0:n] + "..."
	
	return t

## ====================== (!) :


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

# Are those GO IDs in the GO graph?
go_not_in_graph = set(GO2E.keys()) - set(pickle.load(open(IFILE['GO graph'], 'rb')).nodes())
print("Note: {} GO IDs are not in the graph".format(len(go_not_in_graph)))


## DEBUG

#G = pickle.load(open(IFILE['GO graph'], 'rb'))
#print(G['GO:0031625'])
#exit()
#del G


## ===================== WORK :

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

	def subgraph(self, i, depth=8, node_limit=1000) :
		I = [i] # Nodes to be added
		g = nx.MultiDiGraph() # Temporary
		while (I and (depth > 0) and (g.number_of_nodes() < node_limit)) :
			g.add_nodes_from(I)
			E = [(i, j) for i in I for j in sorted(self.G.predecessors(i))]
			E = [(i, j) for (i, j) in E if ("obsolete" not in self.G.nodes[j]['name'])]
			I = [j for (i, j) in E if (j not in g.nodes())]
			g.add_edges_from(E)
			depth -= 1

		g = nx.MultiDiGraph(nx.subgraph(self.G, g.nodes()))
		
		# Remove all metainfo from nodes
		# (it may break the call to graphviz routine)
		for n in g.nodes() :
			info = g.nodes[n]
			for (k, v) in info.items() :
				if (k == 'name') : continue
				g.nodes[n][k] = None
		
		return g


def plot_overview() :
	
	plt.close('all')

	plt.figure(figsize=(math.ceil(6*phi), 6), dpi=150)

	# Plot "almost" all GO terms
	plt.semilogx(
		*zip(*[(n, ci) for (n, CI) in N2CI.items() for ci in CI[0:33]]), 
		'.', color='r', markersize=3
	)
	
	plt.ylim((-1, 1))
	
	plt.xlabel("Size of GO category")
	plt.ylabel("Clustering index")
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['cisz'].format(ext=ext))
	
	plt.close('all')


def plot_tree(root, go_filename) :
	
	plt.close('all')

	print(root, "--", GO2T[root], "({} genes)".format(len(GO2E[root])))
	
	# Ontology subtree
	g = gograph(IFILE['GO graph']).subgraph(root)
	#print(nx.info(g))
	
	#print(root, len(GO2E[root]), GO2E[root])
	
	# For the selection of plt.cm.* colormap, see
	# https://matplotlib.org/examples/color/colormaps_reference.html
	
	node_style = {
		'!' : { 'nodelist' : [], 'node_size' : [], 'node_color' : [] },
		'?' : { 'nodelist' : [], 'node_size' : [], 'node_color' : 'y' },
	}
	
	min_node_size = 20
	
	# Colormap
	cmap = plt.cm.cool
	
	# The ENSG IDs for each GO category that are *also* contained in the root GO category
	GO2e = {
		go : GO2E[root] & GO2E[go]
		for go in g.nodes
	}
	
	for go in g.nodes :
		#ci = GO2CI[go]
		#if ci :
			#node_style['!']['nodelist'].append(go)
			#node_style['!']['node_size'].append(min_node_size + len(GO2E.get(go, {})))
			#node_style['!']['node_color'].append(cmap((ci+1)/2))
		#else :
			#node_style['?']['nodelist'].append(go)
			#node_style['?']['node_size'].append(min_node_size)
		
		q = GO2WQ.get(go, None) # Windowed quantile
		
		if q :
			node_style['!']['nodelist'].append(go)
			node_style['!']['node_size'].append(min_node_size + len(GO2E.get(go, {})))
			node_style['!']['node_color'].append(cmap(q))
		else :
			node_style['?']['nodelist'].append(go)
			node_style['?']['node_size'].append(min_node_size)
		
	
	pos = graphviz_layout(g, prog='dot', root=root)
	# Center
	(cx, cy) = np.mean(list(pos.values()), axis=0)
	pos = { i : (x - cx, y - cy) for (i, (x, y)) in pos.items() }
	# Mirror
	pos = { i : (x, -y) for (i, (x, y)) in pos.items() }
	# Rotate
	a = math.pi * (90/180)
	pos = { i : (np.cos(a) * x - np.sin(a) * y, np.sin(a) * x + np.cos(a) * y) for (i, (x, y)) in pos.items() }
	
	plt.figure(figsize=(10*math.sqrt(2), 10), dpi=300)

	plt.axis('off')
	
	for style in node_style.values() :
		nx.draw_networkx_nodes(g, pos, linewidths=0.1, **style)
	
	nx.draw_networkx_edges(g, pos, width=0.1, arrows=False, edge_color='g')
	
	for (go, (x, y)) in pos.items() :
		
		text = "{} ({})".format(abbr(g.nodes[go]['name'], 45), len(GO2E[go]))
		plt.text(x + 3, y, text, fontsize=4, rotation=0, ha='left', va='bottom')
		
		text = "{}, {}/{} genes in root".format(go, len(GO2e[go]), len(GO2E[go]))
		plt.text(x + 3, y, text, fontsize=1, rotation=0, ha='left', va='top')
	
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
	ax = plt.axes([0.05, 0.05, 0.2, 0.15], facecolor='w')
	
	# Plot "almost" all GO terms
	plt.semilogx(
		*zip(*[(n, ci) for (n, CI) in N2CI.items() for ci in CI[0:33]]), 
		'.', color='y', markersize=1, zorder=-10
	)
	
	plt.xlabel("Size", fontsize=6)
	plt.ylabel("Clustering index", fontsize=6)
	
	ax.tick_params(axis='x', labelsize=6)
	ax.tick_params(axis='y', labelsize=6)
	
	ylim = (-1, +1)   # Possible range of the clustering index
	xlim = plt.xlim() # Range of current GO term sizes
	
	#ax.labelsize = 'small'
	go2size  = dict(zip(node_style['!']['nodelist'], node_style['!']['node_size']))
	go2color = dict(zip(node_style['!']['nodelist'], node_style['!']['node_color']))
	for go in g.nodes :
		#if not (go in GO2E) : continue
		#if not (go in GO2CI) : continue
		if (not GO2E[go]) : continue
		plt.semilogx(len(GO2E[go]), GO2CI[go], 'o', color=go2color[go], markersize=math.sqrt(go2size[go]))
	
	#if (root in GO2CI) :
	plt.semilogx(xlim, (GO2CI[root],) * 2, '-k', linewidth=1)
	
	plt.ylim(ylim)
	plt.xlim(xlim)
	
	# SAVE
	
	for ext in PARAM['ext'] :
		plt.savefig(OFILE['tree'].format(go=go_filename, ext=ext))
	
	plt.close()


def plot_all_trees() :
	
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(plot_tree)(root, str(root).replace(':', '-'))
		for root in sorted(PARAM['GO filter'])
	)
	
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(plot_tree)(root, "Top50_" + x)
		for (root, x) in zip(PARAM['Top50'], ascii_lowercase)
	)
	
	Parallel(n_jobs=PARAM['#proc'])(
		delayed(plot_tree)(root, "Last50_" + x)
		for (root, x) in zip(PARAM['Last50'], ascii_lowercase)
	)


###

def main() :
	plot_all_trees()
	plot_overview()


main()
