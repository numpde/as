
# RA, 2018-01-28

# Run as
#    python3 h2*.py

## ================== IMPORTS :

import os
import sys
import math
import pickle
import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress


## ==================== INPUT :

IFILE = {
	'cossim' : "OUTPUT/g_cossim/UV/TCGA-v-BCXX.pkl",
	
	'clinical' : "OUTPUT/b1_clinical/tcga-clinical.pkl",
}


## =================== OUTPUT :

OFILE = {
	'classified' : "OUTPUT/h2_class/TCGA-BCXX/{method}/{geneset}.pdf",
}


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.2)),
	
	# Classification methods
	'Methods' : {
		'ER status' : {
			'key' : 'er_status_by_ihc',
			'sub' : ["Positive", "Negative", "Indeterminate"],
		},
		
		'HER2 status (FISH)' : {
			'key' : 'her2_fish_status',
			'sub' : ["Positive", "Negative", "Equivocal", "Indeterminate"],
		},
		
		'HER2 status (IHC)' : {
			'key' : 'her2_status_by_ihc',
			'sub' : ["Positive", "Negative", "Equivocal", "Indeterminate"],
		},
	}
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Check if pandas series has unique items
def is_unique(S) : return (S.unique().size == S.size)

def nicer(filename) : 
	return filename.replace(":", "")

def abbr(L) :
	if (type(L) is str) : 
		ab = {
			"Positive" : "(+)",
			"Negative" : "(--)",
			"Equivocal" : "(~)",
			"Indeterminate" : "(?)",
		}
		for (a, b) in ab.items() : L = L.replace(a, b)
		return L
	else :
		return [abbr(i) for i in L]

## ====================== (!) :

pass


## ===================== WORK :

# Compute the similarity grouping by class
def classify(M, method) :
	
	# Clinical classification key in datatable
	clinikey = PARAM['Methods'][method]['key']
	# Possible subtypes
	subtypes = PARAM['Methods'][method]['sub']
	
	# Note:
	# M is a (TCGA bulk sample) x (BCXX single cell sample) similarity matrix
	
	# Partial clinical classification of TCGA samples
	C = pickle.load(open(IFILE['clinical'], 'rb'))['patient_barcode']
	# Select the field of the clinical classification
	C = C.transpose()[clinikey].dropna()
	# Classification subtypes are represented
	assert(set(subtypes) == set(C.unique()))
	
	# Aliquot to class
	def a2c(a) : 
		a = a[0:len(C.index[0])]
		if a in C : return C[a]
		return np.nan
	
	# 
	C = pd.Series(M.index, index=M.index).apply(a2c).dropna()

	# Keep only the TCGA samples that have been classified
	M = M.loc[ list(pd.Series(M.index).isin(C.index)), : ]
	
	# Group into a (TCGA subtype) x (BCXX tumor) matrix
	def group(M) :
		
		# Omit records with missing data
		N = M.dropna().astype(float)
		
		# Group by class
		N = N.groupby(
			# Aliquot barcode to class
			lambda a : C[ a ],
			axis=0
		).mean()
		
		# Group by BCXX tumors
		N = N.groupby(
			# Sample "BCXX[LN][_Re]_YY" to cluster "BCXX"
			lambda c : c[:-3],
			axis=1
		).mean()
		
		# Use the order of the subtypes list
		N = N.reindex(pd.Categorical(N.index, categories=subtypes)).sort_index()
		
		return N
	
	# Group into a (TCGA subtype) x (BCXX tumor) matrix
	N = group(M)

	return N

# Plot a (TCGA subtype) x (BCXX tumor) matrix
def plot(N) :
	
	# Normalize column-wise
	for c in N.columns : N[c] /= N[c].sum()
	
	## Normalize row-wise
	#for i in N.index : N.ix[i] /= N.ix[i].sum()
	
	plt.clf()
	plt.imshow(N.as_matrix(), aspect='auto', cmap=plt.cm.Blues)
	plt.xticks(range(len(N.columns)), list(N.columns), rotation=60)
	plt.yticks(range(len(N.index)), abbr(list(N.index)))
	plt.tick_params(axis='x', which='major', labelsize=6)
	plt.tick_params(axis='y', which='major', labelsize=10)
	plt.colorbar()


def COMPUTE() :
	
	all_cossim = pickle.load(open(IFILE['cossim'], 'rb'))['cossim']
	
	for cossim in Progress()(all_cossim) :
		for method in PARAM['Methods'].keys() :
			meta = cossim['meta']
			setid = meta['id']
			
			#print("Gene set:", setid)
			#print("Method:", method)
		
			N = classify(cossim['M'], method)
			
			plot(N)
			plt.title("{} ({})".format(meta['info'], len(meta['set'])), fontsize=6)
			
			plt.ylabel(method + " normed cos-dist to...")
			
			# Filename for figure
			f = OFILE['classified'].format(method=method, geneset=nicer(setid))
			# Create output directories
			os.makedirs(os.path.dirname(f), exist_ok=True)
			# Save figure
			plt.savefig(f)
	

## ==================== ENTRY :

if (__name__ == "__main__") :
	COMPUTE()

