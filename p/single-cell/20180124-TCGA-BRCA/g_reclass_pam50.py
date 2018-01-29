
# RA, 2018-01-28

# Run as
#    python3 g*.py COMPUTE

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
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
	
	# Gene sets
	'sets' : "OUTPUT/f_gene_subsets/subsets.pkl",
}


## =================== OUTPUT :

OFILE = {
	#'cos-img' : "OUTPUT/g_reclass_pam50/TCGA-v-BCXX.pdf",
	
	'similarity' : "OUTPUT/g_reclass_pam50/UV/TCGA-v-BCXX.pkl",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

PARAM = {
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.5)),
}

PAM50 = ["Normal", "LumA", "LumB", "Her2", "Basal"]

TESTMODE = ("TEST" in sys.argv)


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Log which files are opened
logged_open = open
#def logged_open(filename, mode='r', *argv, **kwargs) :
	#print("({}):\t{}".format(mode, filename))
	#import builtins
	#return builtins.open(filename, mode, *argv, **kwargs)

## ====================== (!) :

# Cosine similarity between two pandas series
def cosine_similarity(a, b) :
	def dot(a, b) : return (a * b).sum()
	norms = math.sqrt(dot(a, a) * dot(b, b))
	return 1 - (dot(a, b) / norms)


## ===================== WORK :


# Compute the similarity along a subset of genes
def compute(genes, meta) :
	
	# BCXX data
	BCXX_DATA = pickle.load(logged_open(IFILE['BCXX'], 'rb'))
	
	# TCGA data
	TCGA_DATA = pickle.load(logged_open(IFILE['TCGA'], 'rb'))

	BCXX = BCXX_DATA['X']
	TCGA = TCGA_DATA['X']
	
	# Keep only the gene symbols that are in both
	(TCGA, BCXX) = TCGA.align(BCXX, join='inner', axis=0)
	
	# Restrict to those genes
	if genes :
		TCGA = TCGA.ix[genes]
		BCXX = BCXX.ix[genes]
	
	# Classification, from
	# http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt
	C = TCGA_DATA['subtype']
	# Narrow down samples
	#C = C[ C['tissue_type'] == "tumor" ]
	# Assume uniqueness in the 'aliquot_barcode' columns
	assert(C['aliquot_barcode'].unique().size == C['aliquot_barcode'].size)
	# Check that the PAM50 types are represented as expected
	assert(set(PAM50) == set(C['PAM50']))
	
	# Expect columns to be unambiguously labelled
	assert(TCGA.columns.unique().size == TCGA.columns.size)
	# Keep only the TCGA samples that have been PAM50-classified
	TCGA = TCGA.loc[ :, TCGA.columns.isin(C['aliquot_barcode']) ]
	assert(set(C['aliquot_barcode']) >= set(TCGA.columns))
	
	# Keep only "nontrivial" samples in the BCXX table
	BCXX = BCXX.loc[ :, BCXX.sum(axis=0) != 0 ]
	
	# Reduce data size for testing
	if TESTMODE :
		TCGA = TCGA[ np.random.permutation(TCGA.columns)[0:9] ]
		BCXX = BCXX[ np.random.permutation(BCXX.columns)[0:7] ]
	
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
	
	# Allocate matrix for (TCGA bulk sample) x (BCXX single cell sample)
	M = pd.DataFrame(index=TCGA.columns, columns=BCXX.columns)
	
	# Compute the similarity matrix
	for a in TCGA.columns :
		for b in BCXX.columns :
			M.loc[a, b] = cosine_similarity(TCGA[a], BCXX[b])
	
	# Group into a (TCGA PAM50 type) x (BCXX tumor) matrix
	N = group(M)
	
	return { 'M': M, 'N' : N, 'meta' : meta }
	
	# How to plot:
	#plt.clf()
	#plt.imshow(N.as_matrix(), aspect='auto', cmap=plt.cm.Blues)
	#plt.xticks(range(len(N.columns)), list(N.columns), rotation=60)
	#plt.yticks(range(len(N.index)), list(N.index))
	#plt.tick_params(axis='x', which='major', labelsize=6)
	#plt.tick_params(axis='y', which='major', labelsize=8)
	#plt.ylabel("Classified as")


def COMPUTE() :
	
	# Load gene subsets of interest
	gene_subsets = pickle.load(logged_open(IFILE['sets'], 'rb'))['S']
	
	if TESTMODE : PARAM['#proc'] = 1

	similarity = Parallel(n_jobs=PARAM['#proc'])(
		delayed(compute)(subset['set'], subset)
		for subset in Progress()(gene_subsets)
	)
	
	pickle.dump(
		{
			'similarity' : similarity,
			'script' : THIS,
		},
		logged_open(OFILE['similarity'], 'wb')
	)

## ==================== ENTRY :

if (__name__ == "__main__") :
	if ("COMPUTE" in sys.argv) : COMPUTE()

