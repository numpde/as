
# RA, 2018-01-30

# Run as
#    python3 i*.py TRAIN
#    python3 i*.py PLOT

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
	# Expression data
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	#'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
	
	# TCGA PAM50 classification
	'TCGA-PAM50' : "OUTPUT/e_prepared/UV/tcga.pkl",
}


## =================== OUTPUT :

OFILE = {
	'model' : "OUTPUT/i1_pam50/UV/trained.pkl",
	
	'tsne-plot' : "OUTPUT/i1_pam50/tsne_{set}.{ext}",
	
	'conf-class' : "OUTPUT/i1_pam50/conf_class_{norm}.{ext}",
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
	
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283537/
	# suggests ["ER", "HER2", "AURKA"], or as HGNC approved symbols:
	'SCMGENE' : ["ESR1", "ERBB2", "AURKA"],
	
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
def is_unique(S) : return (S.unique().size == S.size)


## ====================== (!) :

# Read the PAM50 classification of TCGA samples
def get_tcga_pam50_labels() :
	
	# PAM50
	P = pickle.load(open(IFILE['TCGA-PAM50'], 'rb'))['subtype']
	# Identify record by aliquot barcode
	assert(is_unique(P['aliquot_barcode']))
	# Some patients have multiple aliquots though
	assert(not is_unique(P['patient_barcode']))
	# Use the aliquot barcode as index
	P = P.set_index('aliquot_barcode')
	# Select only the PAM50 class column
	P = P['PAM50']
	# Check that the expected subtypes are represented
	assert(set(P.unique()) == set(PARAM['PAM50-types']))
	# Convert to categorical (useful with P.cat.codes)
	P = pd.Series(pd.Categorical(P.values, categories=PARAM['PAM50-types']), index=P.index)
	
	return P


## ========== WORK (PLOTTING) :

# Show and save the class confusion matrices
# (Normalized over the prediction vector)
def make_confusion_plots(M) :
	
	(X, Y, L, P, I) = (M['X'], M['Y'], M['L'], M['P'], M['I'])
	
	plt.figure(figsize=(10, 5))
	
	# Left and right panels
	for (n, (i, t)) in enumerate([(I, "Training set"), (~I, "Test set")]) :
		(x, y, p) = (X[i], Y[i], P[i])
	
		# http://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html
		# C[i, j] = Number of samples from group i predicted as j
		C = confusion_matrix(np.argmax(y, axis=1), np.argmax(p, axis=1))
		# Normalize the prediction vector for each given class
		C = C / C.sum(axis=1, keepdims=True)
		assert(all(abs(1 - C.sum(axis=1)) <= 1e-10))
		
		plt.subplot(1, 2, 1+n).cla()

		# Note the transpose: display predictions along the y-axis
		plt.imshow(C.T, cmap=plt.cm.Blues, origin='lower', vmin=0, vmax=1)

		plt.xlabel("Reference class")
		plt.xticks(range(len(L)), L)
		
		if (n == 0) :
			plt.ylabel("Classified as...")
			plt.yticks(range(len(L)), L)
		else :
			plt.yticks([])
			
		plt.title(t + " ({} samples)".format(sum(i)))
	
	# "pnorm" means normalization over the prediction vector
	plt.savefig(OFILE['conf-class'].format(norm="pnorm", ext="pdf"))
	plt.close()


# Show and save the t-SNE overview of the reference and predicted classes
def make_tsne_plots(M) :
	
	(X, Y, L, P, I) = (M['X'], M['Y'], M['L'], M['P'], M['I'])
	
	# 2d t-SNE embedding
	from sklearn.manifold import TSNE
	x = TSNE(n_components=2, random_state=1).fit_transform(X).T
	
	# Note: M['x'][d][n] is the n-th run of the d-dim t-SNE embedding
	
	# Convert "probabilities" to class code
	y = np.argmax(Y, axis=1)
	p = np.argmax(P, axis=1)
	
	# Separate scatter plots for the training and the test set
	for (i, s) in [(I, "train"), (~I, "valid")] :
		
		plt.figure(figsize=(8, 5))
		
		# Plot handles for reference and predicted classes
		hy = [None] * len(L)
		hp = [None] * len(L)
		
		# Dummy handles to arrange the legend
		h0 = [ plt.scatter(x[0][0], x[1][0], c="white", s=0) ] * len(L)
		
		## The 'other' samples
		#plt.scatter(x[0][~i], x[1][~i], edgecolors="gray", facecolors="", s=30, marker=10, lw=2)
		#plt.scatter(x[0][~i], x[1][~i], edgecolors="gray", facecolors="", s=30, marker=11, lw=2)
		
		plt.scatter(x[0][(y != p) & i], x[1][(y != p) & i], edgecolors="black", facecolors="", s=70, marker='o', lw=0.5)
		
		# Plot the reference classes of samples
		for (n, _) in enumerate(L) :
			hy[n] = plt.scatter(x[0][(y == n) & i], x[1][(y == n) & i], c=("C{}".format(n)), s=30, marker=10)
		
		# Plot the predicted classes of samples
		for (n, _) in enumerate(L) :
			hp[n] = plt.scatter(x[0][(p == n) & i], x[1][(p == n) & i], c=("C{}".format(n)), s=30, marker=11)
		
		plt.legend(h0 + hy + hp, L + (["Ref"] * len(L)) + (["Est"] * len(L)), ncol=3, prop={'size': 8})
		
		plt.xticks([], [])
		plt.yticks([], [])
		
		plt.savefig(OFILE['tsne-plot'].format(set=s, ext="pdf"))
		plt.close()


def PLOT() :
	M = pickle.load(open(OFILE['model'], 'rb'))
	make_confusion_plots(M)
	make_tsne_plots(M)


## ========== WORK (TRAINING) :

# Train the keras classifier
#    X = samples row-wise (channel first)
#    Y = "sample x class" binary matrix
#    L = labels for the classes
def train(X, Y, L) :
	
	from keras.models import Sequential
	from keras.layers import Dense, Dropout, Activation, BatchNormalization
	from keras        import regularizers
	
	from keras.callbacks import LambdaCallback
	
	from sklearn.model_selection import train_test_split
	
	from keras.regularizers import l2 as L2
	import keras.optimizers
	
	# https://stackoverflow.com/questions/39547279/loading-weights-in-th-format-when-keras-is-set-to-tf-format
	keras_backend.set_image_data_format('channels_first')
	
	# Number of classes
	num_classes = Y.shape[1]
	
	def PassiveInput() : return Dropout(0, input_shape=X.shape[1:])
	def ActiveOutput() : return Dense(Y.shape[1], activation='softmax')
	
	model = Sequential([
		PassiveInput(),
		
		BatchNormalization(),
		
		Dense(8*num_classes, activation='softplus', kernel_regularizer=L2(1e-2)),
		Dropout(0.5),
		Dense(6*num_classes, activation='softplus', kernel_regularizer=L2(1e-2)),
		Dropout(0.5),
		Dense(4*num_classes, activation='softplus', kernel_regularizer=L2(1e-2)),
		Dropout(0.5),
		Dense(2*num_classes, activation='softplus', kernel_regularizer=L2(1e-2)),
		
		ActiveOutput(),
	])
	
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
	
	# Partition into train and validation
	validation_split = 0.2
	np.random.seed(42)
	# Indices of the validation set
	I_valid = (np.random.rand(len(Y)) <= validation_split)
	# Indices of the training set
	I_train = ~I_valid
	
	for _ in range(10) :
		
		# Use a separate small validation split for online inspection
		# Due to the outer loop, this is not a true training/validation split
		model.fit(X[I_train], Y[I_train], epochs=1000, validation_split=0.1)
		
		# Package all data for export
		
		M = {
			'X' : X,
			'Y' : Y,
			'L' : L,
			'I' : I_train,
			'P' : model.predict(X),
			'm' : model,
		}
		
		#plt.ion()
		#plt.show()
		#make_confusion_plots(M)
		#plt.pause(0.1)
	
	
	return M


# Prepare the training data X, Y, L
def get_training_data() :
	
	# For an aliquot_barcode 'a', C[a] is its reference PAM50 type
	C = get_tcga_pam50_labels()
	
	# Load TCGA expression data
	A = pickle.load(open(IFILE['TCGA'], 'rb'))['X']
	
	# Gene set of interest (set to None for all genes)
	genes = None or PARAM['PAM50-genes']
	
	# Has a gene subset been specified?
	if genes :
		# Keep only those genes in the expression table
		A = A.loc[ A.index.isin(genes), : ]
		
		# Requested genes not found in the expression table
		gnitt = (set(genes) - set(A.index))
		if gnitt : print("Genes not in the table:", sorted(gnitt))
	else :
		# The gene subset is "all genes"
		genes = sorted(A.index)
	
	# Index the expression table by aliquot_barcode
	A = A.transpose()
	
	# Keep samples with known class only, aligning the index
	(A, C) = A.align(C, join='inner', axis=0)
	
	# Convert to keras-friendly format
	X = A.astype(float).as_matrix()
	Y = to_categorical(C.cat.codes)
	
	# Class labels
	L = list(C.cat.categories)
	
	assert(X.shape[0] == Y.shape[0])
	assert(X.shape[1] == len(genes))
	assert(Y.shape[1] == len(L))
	
	return (X, Y, L)


def TRAIN() :
	M = train(*get_training_data())
	pickle.dump(M, open(OFILE['model'], 'wb'))



## ==================== ENTRY :

if (__name__ == "__main__") :
	if ("TRAIN" in sys.argv) : TRAIN()
	if ("PLOT"  in sys.argv) : PLOT()
	
	# https://github.com/tensorflow/tensorflow/issues/3388#issuecomment-271107725
	keras_backend.clear_session()
