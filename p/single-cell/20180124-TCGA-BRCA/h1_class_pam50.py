
# RA, 2018-01-28

# Run as
#    python3 h1*.py

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
	
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'clinical' : None,
}


## =================== OUTPUT :

OFILE = {
	'classified' : "OUTPUT/h1_class/TCGA-BCXX/PAM50/{geneset}.pdf",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.2)),
	
	# PAM50 order
	'PAM50' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Check if pandas series has unique items
def is_unique(S) : return (S.unique().size == S.size)

def nicer(filename) : 
	return filename.replace(":", "")

## ====================== (!) :

pass


## ===================== WORK :

# Compute the similarity grouping by class
def classify(M) :
	
	# Note:
	# M is a (TCGA bulk sample) x (BCXX single cell sample) similarity matrix
	
	# Partial classification of TCGA samples, from
	# http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt
	C = pickle.load(open(IFILE['TCGA'], 'rb'))['subtype']
	# These are the PAM50 subtypes
	PAM50 = PARAM['PAM50']
	# Assume uniqueness of 'aliquot_barcode's for unambiguous classification
	assert(is_unique(C['aliquot_barcode']))
	# Expect all the PAM50 subtypes to be represented
	assert(set(PAM50) == set(C['PAM50']))
	
	# Keep only the TCGA samples that have been PAM50-classified
	M = M.loc[ M.index.isin(C['aliquot_barcode']), : ]
	# Each sample is now classified
	assert(set(C['aliquot_barcode']) >= set(M.index))
	
	# Group into a (TCGA PAM50 type) x (BCXX tumor) matrix
	def group(M) :
		
		# Omit records with missing data
		N = M.dropna().astype(float)
		
		# Group by PAM50 type
		N = N.groupby(
			# Aliquot barcode to PAM50 type
			lambda a : C[ C['aliquot_barcode'] == a ]['PAM50'].item(),
			axis=0
		).mean()
		
		# Group by BCXX tumors
		N = N.groupby(
			# Sample "BCXX[LN][_Re]_YY" to cluster "BCXX"
			lambda c : c[:-3],
			axis=1
		).mean()
		
		# Use the order of the PAM50 list
		N = N.reindex(pd.Categorical(N.index, categories=PAM50)).sort_index()
		
		return N
	
	# Group into a (TCGA PAM50 type) x (BCXX tumor) matrix
	N = group(M)

	return N

# Plot a (TCGA type) x (BCXX tumor) matrix
def plot(N) :
	
	N = N.drop(index='Normal')
	
	# Normalize column-wise
	for c in N.columns : N[c] /= N[c].sum()
	
	## Normalize row-wise
	#for i in N.index : N.ix[i] /= N.ix[i].sum()
	
	plt.clf()
	plt.imshow(N.as_matrix(), aspect='auto', cmap=plt.cm.Blues)
	plt.xticks(range(len(N.columns)), list(N.columns), rotation=60)
	plt.yticks(range(len(N.index)), list(N.index))
	plt.tick_params(axis='x', which='major', labelsize=6)
	plt.tick_params(axis='y', which='major', labelsize=8)
	plt.colorbar()
	plt.ylabel("Normed cos-sim to...")


def COMPUTE() :
	
	all_cossim = pickle.load(open(IFILE['cossim'], 'rb'))['cossim']
	
	for cossim in Progress()(all_cossim) :
		setid = cossim['meta']['id']
		#print("Gene set:", setid)
	
		N = classify(cossim['M'])
		
		plot(N)
		plt.title(cossim['meta']['info'], fontsize=6)
		plt.savefig(OFILE['classified'].format(geneset=nicer(setid)))
	

## ==================== ENTRY :

if (__name__ == "__main__") :
	COMPUTE()

