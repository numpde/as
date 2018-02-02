
# RA, 2018-01-30

# Run as
#    python3 i2*.py TRAIN

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

from itertools import chain

from keras.utils import to_categorical
from sklearn.metrics import confusion_matrix

from keras import backend as keras_backend

from utils import make_keras_picklable


## ==================== INPUT :

IFILE = {
	# Expression data
	#'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
	
	# Trained keras model
	'predictor' : "OUTPUT/i1_pam50/UV/trained.pkl",
}


## =================== OUTPUT :

OFILE = {
	'pred-sample' : "OUTPUT/i2_pam50_bcxx/pred_sample.{ext}",
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



## ===================== WORK :

# Plot the predictions of the classifier
def plot(P) :
	from scipy.constants import golden as phi

	# Cell ID to tumor ID
	c2t = (lambda c : c[:-3])
	# Group by tumor
	G2H = dict(P.groupby(by=c2t, axis=1).groups)
	
	# Group names
	G = sorted(G2H.keys())
	# Sample groups
	def key(h) : 
		n = np.argmax(P[h].values)
		return (n, -P[h][n])
		
	S = [sorted(G2H[g], key=key) for g in G]
	
	# 
	R = [len(s) for s in S]
	spacing = 0
	R = list(chain.from_iterable([[spacing, b] for b in R]))[1:]
	R = [sum(R)*0.03] + R
	R = np.cumsum(R)
	R = R / max(R)
	
	fig = plt.figure(figsize=(math.ceil(4*phi), 4), dpi=300)

	AX = [
		fig.add_axes([a, 0, b-a, 1])
		for (a, b) in zip(R[0::2], R[1::2])
	]
	
	origin = "upper"

	for (n, ax) in enumerate(AX) :
		
		ax.imshow(
			P[S[n]].as_matrix(), 
			aspect='auto', origin=origin, cmap=plt.cm.Blues,
			vmin=-0.1, vmax=1
		)
		
		ax.set_xticks([])
		ax.set_yticks([])
		
		# Group label
		ax.text(
			0.95, 0.005, G[n], 
			size=6,
			transform=ax.transAxes, rotation=90, 
			ha='right', va='bottom',
		)
	
	ax = fig.add_axes([0, 0, R[0], 1])
	L = list(P.index)
	if (origin == "upper") : L = list(reversed(L))
	for (y, t) in zip(np.linspace(0, 1, 1 + 2*len(L))[1::2], L) :
		#ax.plot(0.1, y, '.')
		#ax.set_xlim(0, 1)
		#ax.set_ylim(0, 1)
		ax.text(
			0.5, y, t, 
			size=8,
			transform=ax.transAxes, rotation=90, 
			ha='center', va='center',
		)
		
	
	#plt.imshow(P, cmap=plt.cm.Blues, origin='lower', vmin=0, vmax=1, aspect='auto')
	#plt.ylabel("Classified as...")
	#plt.yticks(range(len(P.index)), P.index)
	##plt.xticks(range(len(P.columns)), P.columns)
	
	plt.savefig(OFILE['pred-sample'].format(ext="pdf"))
	
	#plt.show()


def predict_bcxx() :
	# BCXX expression as pandas table
	B = pickle.load(open(IFILE['BCXX'], 'rb'))['X']
	
	## Cell ID to tumor ID
	#c2t = (lambda c : c[:-3])
	### Classify the average expression in each tumor
	##B = B.groupby(by=c2t, axis=1).mean()
	
	# Drop invalid data
	B = B.dropna()
	
	# Trained PAM50 predictor bundle
	M = pickle.load(open(IFILE['predictor'], 'rb'))
	
	# keras model
	m = M['m']
	
	# Select genes = features
	B = B.loc[ M['F'], : ]
	
	# Normalize sample-wise
	for c in B.columns : B[c] /= B[c].sum()
	
	# Drop missing/invalid data
	B = B.dropna()

	# Class labels
	L = M['L']
	
	# Assemble sample classifiction into a dataframe
	P = pd.DataFrame(m.predict(B.T).T, index=L, columns=B.columns)
	
	return P


## ==================== ENTRY :

if (__name__ == "__main__") :
	plot(predict_bcxx())
	
	# https://github.com/tensorflow/tensorflow/issues/3388#issuecomment-271107725
	keras_backend.clear_session()
