
# RA, 2018-01-28

# Run as
#    python3 g*.py

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
	'cossim' : "OUTPUT/g_cossim/UV/TCGA-v-BCXX.pkl",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.5)),
	
	# Record just in case
	'testmode' : TESTMODE,
	
	# Do gene-wise zscore transform?
	'zscore' : False,
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Check if pandas series has unique items
def is_unique(S) : return (S.unique().size == S.size)

# zscore transform of a pandas series
def zscore(S) : return (S - S.mean()) / S.std()


## ====================== (!) :

# Cosine similarity between two pandas dataframes (column-wise)
def cosine_similarity(A, B) :
	
	# Drop zero columns
	A = A.loc[ :, (A != 0).any(axis=0) ]
	B = B.loc[ :, (B != 0).any(axis=0) ]
	
	# Convert to numpy arrays
	a = A.astype(float).as_matrix()
	b = B.astype(float).as_matrix()
	
	# All-x-all dot products
	d = np.matmul(a.T, b)
	
	# Norm-inverse of columns
	na = np.diag(1 / np.linalg.norm(a, axis=0))
	nb = np.diag(1 / np.linalg.norm(b, axis=0))
	
	# Cosines
	c = na.dot(d).dot(nb)
	
	# Cosine similarities
	s = 1 - c
	del c
	
	# Check expected dimensions
	assert(s.shape == (len(A.columns), len(B.columns)))
	
	# Return as dataframe
	return pd.DataFrame(s, index=A.columns, columns=B.columns)


## ===================== WORK :

# Compute the similarity along a subset of genes
def compute(genes, meta) :
	
	# BCXX data
	BCXX_DATA = pickle.load(open(IFILE['BCXX'], 'rb'))
	
	# TCGA data
	TCGA_DATA = pickle.load(open(IFILE['TCGA'], 'rb'))

	BCXX = BCXX_DATA['X']
	TCGA = TCGA_DATA['X']
	
	# Restrict to those genes
	if genes :
		TCGA = TCGA.ix[genes]
		BCXX = BCXX.ix[genes]
	
	if PARAM['zscore'] :
		TCGA = TCGA.apply(zscore, axis=1)
		BCXX = BCXX.apply(zscore, axis=1)
	
	TCGA = TCGA.dropna()
	BCXX = BCXX.dropna()
	
	# Keep only the genes that are in both tables
	# Note: the tables are indexed by the gene symbol
	(TCGA, BCXX) = TCGA.align(BCXX, join='inner', axis=0)
	
	# Similarity matrix for (TCGA bulk sample) x (BCXX single cell sample)
	M = cosine_similarity(TCGA, BCXX)
	
	return { 'M' : M, 'meta' : meta }


def COMPUTE() :
	
	# Load gene subsets of interest
	gene_subsets = pickle.load(open(IFILE['sets'], 'rb'))['S']

	cossim_by_geneset = Parallel(n_jobs=PARAM['#proc'])(
		delayed(compute)(geneset['set'], geneset)
		for geneset in Progress()(gene_subsets)
	)
	
	pickle.dump(
		{
			'cossim' : cossim_by_geneset,
			'script' : THIS,
			'param'  : PARAM,
		},
		open(OFILE['cossim'], 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	COMPUTE()

