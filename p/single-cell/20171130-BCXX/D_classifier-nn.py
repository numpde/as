
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
	
	'tsne runs' : "OUTPUT/C_goordinates/tsne_runs.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'classified' : "OUTPUT/D_classifier-nn/classified_epochs={epochs}.pkl",
	
	'samples'    : "OUTPUT/D_classifier-nn/samples.{ext}",
	'confusion'  : "OUTPUT/D_classifier-nn/confusion.{ext}",
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

# This script
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

## ====================== (!) :

pass

## ===================== DATA :

#[ BC DATA ]#

# Load the BC data
BC_data = pickle.load(open(IFILE['BC data'], 'rb'))

# Expression matrix
BC_X = BC_data['X']

# Rearrange data axes
(axis_smpl, axis_gene) = (0, 1)
BC_X = np.moveaxis(BC_X, BC_data['axis_gene'], axis_gene)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (BC_X.shape[axis_smpl], BC_X.shape[axis_gene])

# ENSG IDs
BC_E = BC_data['gene_id']
assert(len(BC_E) == BC_X.shape[axis_gene]), "Inconsistent gene info"

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


## The Gene Ontology graph
#GO_graph = pickle.load(open(IFILE['GO graph'], 'rb'))

## Are those GO IDs in the GO graph?
#go_not_in_graph = set(GO2E.keys()) - set(GO_graph.nodes())
#print("Note: {} GO IDs are not in the graph".format(len(go_not_in_graph)))


## ===================== WORK :

#[ ]#


def test() :
	go = "GO:0006260" # DNA replication
	X = np.take(BC_X, GO2I[go], axis=axis_gene)
	
	#rec = [rec for rec in pickle.load(open(IFILE['tsne runs'], 'rb'))['runs'] if ((rec['N'] == 200) and (rec['run'] == 1))][0]
	#X = np.moveaxis(rec['Y'], rec['axis_smpl'], axis_smpl)
	#assert(X.shape[axis_smpl] == n_samples)
	
	assert(axis_smpl == 0)
	
	## Normalize gene-wise
	#assert(axis_gene == 1)
	#for n in range(X.shape[1]) : X[:, n] /= np.sum(X[:, n])
	
	## Normalize sample-wise
	#assert(axis_smpl == 0)
	#for c in range(X.shape[0]) : X[c, :] /= np.sum(X[c, :])
	
	num_classes = len(G2S)
	Y = np.zeros((X.shape[axis_smpl], num_classes))
	
	# Reference classes
	for (n, s) in enumerate(S) : Y[s, n] = 1
	
	#plt.imshow(Y); plt.show(); exit()

	from keras.models import Sequential
	from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten, Activation
	from keras        import regularizers
	
	from keras.callbacks import LambdaCallback

	# https://stackoverflow.com/questions/39547279/loading-weights-in-th-format-when-keras-is-set-to-tf-format
	from keras import backend as keras_backend
	keras_backend.set_image_dim_ordering('th')
	
	import keras.regularizers
	import keras.optimizers
	
	model = Sequential([
		Dropout(rate=0.6, input_shape=(X.shape[axis_gene],)),
		Dense(num_classes, kernel_regularizer=keras.regularizers.l1_l2(1, 1)),
		Activation('softmax'),
	])
	
	opt = 'adam'
	#opt = keras.optimizers.SGD(lr=0.0001, decay=1e-6, momentum=0.9, nesterov=True)
	model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
	#model.compile(loss='mean_absolute_error', optimizer=opt, metrics=['accuracy'])
	
	#figs = { 1 : plt.figure(), 2 : plt.figure() }
	
	def plot_model_prediction() :
		
		Yp = model.predict(X)
		
		def sort_by_prediction(Yp) :
			for (n, s) in enumerate(S) :
				Yp[s, :] = Yp[sorted(s, key=(lambda c : -Yp[c, n])), :]
			return Yp
	
		plt.ion()
		
		plt.figure(1)
		
		f = plt
		
		f.imshow(sort_by_prediction(Yp).T, aspect='auto')
		f.xlabel("Samples")
		f.ylabel("Classified as")
		f.gca().xaxis.set_label_position('top') 
		f.xticks([], [])
		f.yticks(range(num_classes), sorted(G2S.keys()), rotation=45)
		f.tick_params(axis='both', which='major', labelsize=8)
		for ext in PARAM['ext'] : f.savefig(OFILE['samples'].format(ext=ext))
		#f.show()
		
		# Confusion matrix
		C = np.vstack(
			np.mean(Yp[s, :], axis=0)
			for (g, s) in sorted(G2S.items())
		)
		
		plt.figure(2)
		
		f.imshow(C.T, vmin=0, vmax=1)
		f.xlabel("Class")
		f.ylabel("Classified as")
		f.gca().xaxis.set_label_position('top') 
		f.xticks(range(num_classes), sorted(G2S.keys()), rotation=45)
		f.yticks(range(num_classes), sorted(G2S.keys()), rotation=45)
		f.tick_params(axis='both', which='major', labelsize=8)
		for ext in PARAM['ext'] : f.savefig(OFILE['confusion'].format(ext=ext))
		#f.show()
		
		plt.pause(0.1)
	
	
	def save_model_prediction(epochs) :
		pickle.dump(
			{ 'Yp' : model.predict(X), 'script' : THIS }, 
			open(OFILE['classified'].format(epochs=epochs), 'wb')
		)
	
	def every_epoch(epoch, logs) :
		
		epochs = epoch + 1
		
		if (0 == (epochs % 1000)) :
			plot_model_prediction()
		
		if (0 == (epochs % 1000)) :
			save_model_prediction(epochs)
	
	model.fit(X, Y, epochs=9000, callbacks=[LambdaCallback(on_epoch_end=every_epoch)])
	
	##
	
	

	input()


###

if (__name__ == "__main__") :
	#plot_CI_vs_N()
	test()
