
# RA, 2018-01-05

# IN PROGRESS

# The GO category
# 	http://amigo.geneontology.org/amigo/term/GO:0032201
# 	telomere maintenance via semi-conservative replication
# shows a relatively low clustering index and
# a high degree of intra-tumoral heterogeneity
#
# Heatmap:
#	~/svn/as/p/single-cell/20171130-BCXX/OUTPUT/7_heatmaps/heatmap_GO-0032201.png
# GO context:
#	~/svn/as/p/single-cell/20171130-BCXX/OUTPUT/9_tree/tree_Top50_c.png
#
# Here, try to identify the clusters in this GO category

# Ref:
# https://bib.dbvis.de/uploadedFiles/MatrixReorderingSTAR.pdf
# https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
# https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html

## ================ IMPORTS 1 :

import matplotlib
#
# https://stackoverflow.com/questions/15417586/python-matlplotlib-add-hyperlink-to-text/30099364#30099364
matplotlib.use("pgf")
pgf_with_custom_preamble = {
	"text.usetex": True,
	"pgf.preamble": [ r"\usepackage[hidelinks]{hyperref}" ],
}
matplotlib.rcParams.update(pgf_with_custom_preamble)

## ================ IMPORTS 2 :

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from sklearn.decomposition import PCA

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
from scipy.cluster.hierarchy import linkage
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.leaves_list.html
from scipy.cluster.hierarchy import leaves_list

## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'GO=>Info' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'GO=>CI'   : "OUTPUT/0_go2ci/UV/go2ci.pkl",
	'gene-ks'  : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'clustering1' : "OUTPUT/A_followup/clustering_{scope}_{go}.pdf",
	'cross-go'    : "OUTPUT/A_followup/cross-go_{go}_{group}.pdf",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	'GO (cluster)' : {
		'GO:0032201', # telomere maintenance via semi-conservative replication
		'GO:0010833', # telomere maintenance via telomere lengthening
		'GO:0005333', # norepinephrine transmembrane transporter activity
		
		'GO:0030284', # estrogen receptor activity
		'GO:0030520', # intracellular estrogen receptor signaling pathway
		'GO:0038128', # ERBB2 signaling pathway (ERBB2 = ENSG00000141736)
	},
	
	
	'GO (cross)' : [
		[
			# Compare expression of
			'GO:0032201', # telomere maintenance via semi-conservative replication
			# against
			'GO:0010833', # telomere maintenance via telomere lengthening
		],
		
		[
			# Compare expression of
			"GO:0035004", # phosphatidylinositol 3-kinase activity
			# against
			"GO:0051019", # mitogen-activated protein kinase binding
		],	
	]
}


## ================== HELPERS :

def safe(go) :
	return str(go).replace(':', '-')

def ensg_url(e) :
	return "\\href{https://asia.ensembl.org/Homo_sapiens/Gene/Summary?g=" + e + "}{" + e + "}"

## ===================== DATA :


#[ LOAD BC DATASET ]#

# Load the BC data
BC_data = pickle.load(open(IFILE['BC data'], 'rb'))

# Expression matrix
BC_X = BC_data['X']

# Labels for axis/dimension of BC data
(axis_smpl, axis_gene) = (BC_data['axis_smpl'], BC_data['axis_gene'])

# ENSG IDs
BC_E = BC_data['gene_id']
assert(len(BC_E) == BC_X.shape[axis_gene])

# E2I : ENSG ID --> Index
E2I = dict(zip(BC_E, range(len(BC_E))))

# Clusters/groups
G2S = { 
	g : tuple(s for (s, h) in SH)
	for (g, SH) in BC_data['B2SH'].items() 
}
S = list(G2S.values())

# I2G : Sample index --> Group number
I2G = { i : n  for (n, s) in enumerate(sorted(S))  for i in s }


#[ LOAD GO TERMS ]#

# GO2E : GO ID --> [ENSG IDs]
GO2E = {
	go_E[0] : set(go_E[1:])
	for go_E in [
		L.rstrip().split('\t') 
		for L in open(IFILE['GO=>ENSG'], 'r')
	]
}


#[ LOAD KS DIFFERENTIAL EXPRESSION DATA ]#

KS_data = pickle.load(open(IFILE['gene-ks'], 'rb'))
#print(KS.keys())

KS_OP = 1
KS_meta = KS_data['KS_meta'][KS_OP]
assert(KS_meta in ['max', 'mean', 'median', 'min', 'std'])

# Differential expression by ENSG
E2DE = { e : ks[KS_OP] for (e, ks) in KS_data['E2KS'].items() }


## ===================== WORK :

def plot_clustering() :
	
	# GO term of interest
	for go in PARAM['GO (cluster)'] :
		
		for scope in ['global', 'local'] :

			# Select the genes from the GO category of interest
			# Sort them by differential expression
			E = list(sorted(GO2E[go], key=(lambda e : -E2DE[e])))

			# BC data subset for those genes
			X = np.take(BC_X, [E2I[e] for e in E], axis=axis_gene)

			# Z-transform gene-wise
			Z = stats.mstats.zscore(X, axis=axis_smpl)

			# Number of samples / genes in the reduced expression matrix
			(n_samples, n_genes) = (Z.shape[axis_smpl], Z.shape[axis_gene])

			print("{}: {} genes. Clustering scope: {}.".format(go, n_genes, scope))

			#[ PCA ]#

			# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html

			#pca = PCA(n_components=(n_genes-1))
			#pca.fit(Y)
			#Z = pca.transform(Y)
			#assert(Z.shape[0] == n_samples)

			##plt.plot(Z[:, 1], Z[:, 2], 'x')
			##plt.show()

			#plt.imshow(Z)
			#plt.show()


			#[ HIERARCHICAL CLUSTERING ]#

			# Compute a distance matrix as (1 - cos(angle))
			def cos_dist_mat(X, axis) :
				# Covariance & norm products
				C = np.tensordot(X, X, axes=([axis], [axis]))
				V = np.sqrt(np.outer(np.diag(C), np.diag(C)))
				V[V == 0] = 1
				D = np.max(1 - (C / V), 0)
				return D

			def cos_dist(x, y) :
				return 1 - (np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y))

			assert(axis_smpl == 0), "The implementation is not general"
			
			if (scope == 'global') :
				L = linkage(Z, method='complete') #, metric=cos_dist)
				#
				#Y = np.zeros(X.shape)
				#Y[X > 1] = 1
				#L = linkage(Y, method='complete')
				#

				#plt.imshow(cos_dist(Y, axis_gene))
				#plt.show()

				J = list(leaves_list(L))
				del L
			
			if (scope == 'local') :
				J = []
				
				for s in sorted(S) :
					z = np.take(Z, s, axis=axis_smpl)
					L = linkage(z, method='complete') #, metric=cos_dist)
					#
					#Y = np.zeros(X.shape)
					#Y[X > 1] = 1
					#L = linkage(Y, method='complete')
					#

					#plt.imshow(cos_dist(Y, axis_gene))
					#plt.show()

					J.extend(list(s[i] for i in leaves_list(L)))
					del L
			
			plt.close('all')
			plt.clf()
			
			(f, AX) = plt.subplots(4)
			#
			for ax in AX : 
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				ax.xaxis.set_ticks_position('none')
				ax.yaxis.set_ticks_position('none')
				
				# https://stackoverflow.com/questions/1639463/matplotlib-border-width
				[i.set_linewidth(0.1) for i in ax.spines.values()]
			#
			#
			ax = AX[0]
			XX = np.zeros(X.shape)
			XX[X >= 1] = np.log(X[X >= 1])
			ax.imshow(XX.transpose(), aspect='auto')
			ax.yaxis.set_ticks(range(0, n_genes))
			ax.yaxis.set_ticklabels(map(ensg_url, E))
			[tick.label.set_fontsize(1) for tick in ax.yaxis.get_major_ticks()]
			#list(range(1, 1 + n_genes)), list(E))
			#
			ax = AX[1]
			ax.imshow(np.asarray([(I2G[i], I2G[i]) for i in range(n_samples)]).transpose(), aspect='auto')
			#
			#
			ax = AX[2]
			XX = np.zeros(X.shape)
			XX[X >= 1] = np.log(X[X >= 1])
			XX = XX[J, :]
			ax.imshow(XX.transpose(), aspect='auto')
			ax.yaxis.set_ticks(range(0, n_genes))
			ax.yaxis.set_ticklabels(map(ensg_url, E))
			[tick.label.set_fontsize(1) for tick in ax.yaxis.get_major_ticks()]
			#
			ax = AX[3]
			ax.imshow(np.asarray([(I2G[j], I2G[j]) for j in J]).transpose(), aspect='auto')
			#
			#
			plt.savefig(OFILE['clustering1'].format(scope=scope[0:3], go=safe(go)))


def plot_crossgo() :
	
	for GO in PARAM['GO (cross)'] :
		
		assert(type(GO) is list)
		assert(len(GO) == 2)
		
		print(GO[0], "vs", GO[1])

		# Select the genes from the GO category of interest
		# Sort them by differential expression
		E = {
			go : list(sorted(GO2E[go], key=(lambda e : -E2DE[e])))
			for go in GO
		}
		
		# The two sets of genes are almost disjoint!
		assert(len(set(E[GO[0]]) & set(E[GO[1]])) <= 1)

		# BC data subset for those genes
		X = {
			go : np.take(BC_X, [E2I[e] for e in E[go]], axis=axis_gene)
			for go in GO
		}

		for g0 in sorted(G2S.keys()) :
			
			plt.close('all')
			plt.clf()
			
			# Handles for the legend
			H1 = [] 
			H2 = []
			
			for (g, s) in sorted(G2S.items()) :
				
				def nz(x) :
					#if any(x == 0) :
					x[x == 0] = np.min(x[x > 0] / 10)
					return x
				
				x = [
					nz(np.mean(np.take(X[go], s, axis=axis_smpl), axis=axis_gene))
					for go in GO
				]
				
				if (g == g0) : 
					c = 'b'
					z = 0
				else :
					c = 'y'
					z = -10

				H1.append(plt.loglog(*x, 'o', color=c, zorder=z)[0])
				
				if (g != g0) : continue
			
				def prox(i) : return np.log(x[0][i] / x[1][i])
				
				a = max(range(len(s)), key=prox)
				b = min(range(len(s)), key=prox)
				H2.append(plt.loglog((x[0][a], x[0][b]), (x[1][a], x[1][b]), ':b')[0])
				
			xlim = plt.xlim()
			ylim = plt.ylim()
			ones = (1, 1)
			
			plt.loglog(xlim, xlim, '--k', linewidth=0.5)
			plt.loglog(ones, ylim, '--k', linewidth=0.5)
			plt.loglog(xlim, ones, '--k', linewidth=0.5)
			
			L1 = sorted(G2S.keys())
			L2 = ["distant"]
			
			H1.extend(H2)
			L1.extend(L2)
			
			plt.legend(H1, L1, prop={'size': 6}, loc='lower left')
			
			
			plt.xlabel("Mean expression in " + GO[0])
			plt.ylabel("Mean expression in " + GO[1])
			
			goxgo = "{}x{}".format(safe(GO[0]), safe(GO[1]))
			plt.savefig(OFILE['cross-go'].format(go=goxgo, group=g0))

# # #

any_work_done = False

if ('1' in sys.argv) : 
	plot_clustering()
	any_work_done = True

if ('2' in sys.argv) : 
	plot_crossgo()
	any_work_done = True

if not any_work_done :
	print("Specify switch '1' (for clustering) and/or '2' (for cross-go)")
