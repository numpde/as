
# RA, 2018-01-30

# Run as
#    python3 i3*.py

# Requires the module
#    utils/make_keras_picklable.py


## ================== IMPORTS :

import os
import sys
import math
import pickle
import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from keras.utils import to_categorical
from sklearn.metrics import confusion_matrix

from keras import backend as keras_backend

from utils import make_keras_picklable


## ==================== INPUT :

IFILE = {
	# METABRIC expression data (log-intensity)
	'BRIC' : "OUTPUT/e_prepared/UV/bric.pkl",
	# TCGA expression data (FPKM)
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	
	# TCGA-trained predictor
	'model' : "OUTPUT/i1_pam50/UV/trained.pkl",
}


## =================== OUTPUT :

OFILE = {
	#'model' : "OUTPUT/i1_pam50/UV/trained.pkl",
	
	#'tsne-plot' : "OUTPUT/i1_pam50/tsne_{set}.{ext}",
	
	'conf-class' : "OUTPUT/i3_pam50_bric/conf_class_{norm}.{ext}",
	
	#'model-plot' : "OUTPUT/i1_pam50/model.{ext}",
	
	#'expression-heatmap' : "OUTPUT/i1_pam50/expression-heatmap/{pam50}.{ext}",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	## Number of parallel computing processes
	#'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.2)),
	
	# Record just in case
	'testmode' : TESTMODE,
	
	# Source: https://www.biostars.org/p/77590/
	# Replaced ORC6L by ORC6 (http://www.genecards.org/cgi-bin/carddisp.pl?gene=ORC6)
	'PAM50-genes' : ("ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 NDC80 NUF2 ORC6 PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T").split(),
	
	'PAM50-types' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
	
	## Extension for plots
	#'ext' : { 'pdf', 'png' },
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Check if pandas series has unique items
def is_unique(S) : return (len(set(S)) == len(S))

def restrict_to_pam50(X) :
	return X.ix[ X.index.isin(PARAM['PAM50-genes']) ]

## ====================== (!) :


# Normalization as in the training
def normalize(A) :
	
	# Normalize columns-wise
	for c in A.columns : A[c] /= (sum(A[c]) or 1)
	
	# Log-transform
	A = np.log(1 + 10 * A)
	
	return A


def predict(M, B) :
	
	# Classify samples and assemble into a dataframe
	# Normalization should be done only once!
	P = pd.DataFrame(M['m'].predict(normalize(B).as_matrix().T).T, index=M['L'], columns=B.columns)

	return P


# Adapt row-wise mean and variance of B to that of A
def mean_std_normalize_row(B, A) :
	
	def ascore(r) :
		s = A.ix[r.name]
		return np.mean(s) + np.std(s) * ((r - np.mean(r)) / np.std(r))
	
	return B.apply(ascore, axis=1)


# Column-wise quantile normalization of B using A
# https://en.wikipedia.org/wiki/Quantile_normalization
def quantile_normalize_col(B, A) :
	
	# A and B are of layout [features] x [samples]
	
	# Get A's rows for normalization
	(B, A) = B.align(A, join='left', axis=0)
	
	# v[r] is the typical value of rank r
	v = np.mean(np.sort(A.as_matrix(), axis=0), axis=1)
	
	for c in B : B[c][np.argsort(B[c].values)] = v
	
	assert(not B.isnull().values.any()), "Behavior with NaN is undefined"
	
	return B


## ========== WORK (PLOTTING) :

def get_bric_data() :
	
	BRIC = pickle.load(open(IFILE['BRIC'], 'rb'))
	TCGA = pickle.load(open(IFILE['TCGA'], 'rb'))
	
	# Predictor
	M = pickle.load(open(IFILE['model'], 'rb'))
	
	X = BRIC['X'] # Expression matrix
	C = BRIC['C'] # Classification
	
	#print(set(X.index[ X.index.isin(["BAG1", "BAG-1", "HAP", "RAP46"]) ]))
	#print(set(X.index[ X.index.isin(["GPR160", "GPCR150", "HGPCR1", "GPCR1", "23693"]) ]))
	#print(set(X.index[ X.index.isin(["MIA", "CD-RAP"]) ]))
	#print(set(X.index[ X.index.isin(["TMEM45B", "Transmembrane Protein 45B"]) ]))
	
	X = restrict_to_pam50(X)
	
	G = sorted(set(PARAM['PAM50-genes']) - set(X.index))
	print("PAM50   genes not in METABRIC:", (", ".join(G)))
	
	G = sorted(set(M['F']) - set(X.index))
	print("Trained genes not in METABRIC:", (", ".join(G)))
	
	# log-intensity to (sort of) relative counts
	X = np.exp(np.log(2) * X)
	
	
	# Adapt gene-wise mean/variance to TCGA
	X = mean_std_normalize_row(X, TCGA['X'])
	X[ X < 0 ] = 0
	
	# Column-wise quantile normalization
	#X = quantile_normalize_col(X, pickle.load(open(IFILE['TCGA'], 'rb'))['X'])
	
	
	# Populate with any missing genes according to the training phase
	X = X.append( pd.DataFrame(0, index=G, columns=X.columns) ).sort_index()
	
	# Reference class
	Y = pd.DataFrame(list(C['CLAUDIN_SUBTYPE']), index=C['PATIENT_ID'], columns=['PAM50'])
	Y = Y[ Y['PAM50'].isin(PARAM['PAM50-types']) ]
	
	# Align by sample (Y.index and X.columns)
	Y = Y[ Y.index.isin(X.columns) ]
	X = X[ Y.index ]
	
	# Convert to categorical datatype
	Y['PAM50'] = pd.Series(pd.Categorical(Y['PAM50'].values, categories=PARAM['PAM50-types']), index=Y.index)
	
	Y = to_categorical(Y['PAM50'].cat.codes, len(Y['PAM50'].cat.categories))
	
	return (X, Y)


# Show and save the class confusion matrices
# (Normalized over the prediction vector)
def plot_confusion_matrices() :
	
	# Predictor
	M = pickle.load(open(IFILE['model'], 'rb'))
	
	(X, Y) = get_bric_data()
	
	P = predict(M, X).as_matrix().T
	L = PARAM['PAM50-types']
	
	plt.figure(figsize=(5, 5), dpi=100)
	
	# Convert "probabilities" to class code
	y = np.argmax(Y, axis=1)
	p = np.argmax(P, axis=1)
	
	# Left and right panels
	
	t = "METABRIC"
	
	# http://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html
	# C[i, j] = Number of samples from group i predicted as j
	C = confusion_matrix(y, p)
	# Normalize the prediction vector for each given class
	C = C / C.sum(axis=1, keepdims=True)
	assert(all(abs(1 - C.sum(axis=1)) <= 1e-10))
	
	#plt.subplot(1, 2, 1+n).cla()

	# Note the transpose: display predictions along the y-axis
	im = plt.imshow(C.T, cmap=plt.cm.Blues, origin='lower', vmin=0, vmax=1)
	ax = plt.gca()
	
	ax.tick_params(axis='both', which='both', length=0)
	
	plt.xticks(range(len(L)), L)
	plt.xlabel("Reference class")
	plt.yticks(range(len(L)), L)
	plt.ylabel("Classified as...")
	
	plt.title(t + " ({}/{} samples correct)".format(sum(y == p), len(y)))
	
	plt.tight_layout()
	
	(y0, y1) = (ax.get_position().y0, ax.get_position().y1)
	
	# Colorbar
	
	# https://stackoverflow.com/questions/13784201/matplotlib-2-subplots-1-colorbar
	cax = plt.gcf().add_axes([0.98, y0, 0.01, y1-y0])
	plt.gcf().subplots_adjust(right=0.9)
	cb = plt.colorbar(im, cax=cax, ticks=[0, 0.5, 1])
	cax.tick_params(labelsize=5) 
	#
	cb.set_clim(0, 1)
	cax.set_yticklabels(["0%", "proportion", "100%"], va='center', rotation=90)
	cax.yaxis.set_ticks_position('left')
	cax.tick_params(axis='both', which='both', length=0)
	
	# Save figure to disk
	
	# "pnorm" means normalization over the prediction vector
	plt.savefig(OFILE['conf-class'].format(norm="pnorm", ext="pdf"))
	plt.close()



def PLOT() :
	plot_confusion_matrices()



## ==================== ENTRY :

if (__name__ == "__main__") :
	PLOT()
	
	# https://github.com/tensorflow/tensorflow/issues/3388#issuecomment-271107725
	keras_backend.clear_session()
