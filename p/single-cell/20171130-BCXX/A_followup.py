
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

## ================== IMPORTS :

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from sklearn.decomposition import PCA


## ==================== INPUT :

IFILE = {
	'BC data'  : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
	'GO=>ENSG' : "OUTPUT/0_e2go/go2e.txt",
	'GO=>Info' : "OUTPUT/0_go-graph/UV/go-graph.pkl",
	'GO=>CI'   : "OUTPUT/0_go2ci/UV/go2ci.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)


## ==================== PARAM :

PARAM = {
	'GO' : 'GO:0032201', # telomere maintenance via semi-conservative replication
}


## ===================== DATA :


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
assert(len(BC_E) == n_genes)

# E2I : ENSG ID --> Index
E2I = dict(zip(BC_E, range(len(BC_E))))

#[ LOAD GO TERMS ]#

# GO2E : GO ID --> [ENSG IDs]
GO2E = {
	go_E[0] : set(go_E[1:])
	for go_E in [
		L.rstrip().split('\t') 
		for L in open(IFILE['GO=>ENSG'], 'r')
	]
}


## ===================== WORK :

# Z-transform gene-wise
X = stats.mstats.zscore(X, axis=axis_smpl)

# GO term of interest
go = PARAM['GO']

# Select the genes from the GO category of interest
Y = np.take(X, [E2I[e] for e in GO2E[go]], axis=axis_gene)
del X

# Number of samples / genes in the reduced expression matrix
(n_samples, n_genes) = (Y.shape[axis_smpl], Y.shape[axis_gene])

print("{} genes left.".format(n_genes))


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
def cos_dist(X, axis) :
	# Covariance & norm products
	C = np.tensordot(X, X, axes=([axis], [axis]))
	V = np.sqrt(np.outer(np.diag(C), np.diag(C)))
	V[V == 0] = 1
	D = 1 - (C / V)
	return D

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
from scipy.cluster.hierarchy import linkage
assert(axis_smpl == 0)
Z = linkage(Y, method='complete', metric=(lambda x, y : cos_dist(np.vstack([x, y]), axis=1)[0, 1]))
print(Z)

#plt.imshow(cos_dist(Y, axis_gene))
#plt.show()

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.leaves_list.html
from scipy.cluster.hierarchy import leaves_list

plt.imshow(cos_dist(Y[np.asarray(leaves_list(Z)), :], axis_gene))
plt.show()

