
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
# https://en.wikipedia.org/wiki/Cosine_similarity
def cosine_similarity(A, B) :

	# Nonzero columns
	jA = (A != 0).any(axis=0)
	jB = (B != 0).any(axis=0)

	# Keep only nonzero columns
	A = A.loc[ :, jA ]
	B = B.loc[ :, jB ]
	
	# Convert to numpy arrays
	a = A.astype(float).values
	b = B.astype(float).values

	# All-x-all dot products
	d = np.matmul(a.T, b)
	
	# Norm-inverse of columns
	na = np.diag(1 / np.linalg.norm(a, axis=0))
	nb = np.diag(1 / np.linalg.norm(b, axis=0))
	
	# Cosines
	c = na.dot(d).dot(nb)
	
	# Check expected dimensions
	assert(c.shape == (len(A.columns), len(B.columns)))
	
	# Return as dataframe
	return pd.DataFrame(c, index=A.columns, columns=B.columns)


## ===================== WORK :

# Compute the similarity along a subset of genes
def compute(genes, meta) :
	
	# BCXX data
	BCXX_DATA = pickle.load(open(IFILE['BCXX'], 'rb'))
	
	# TCGA data
	TCGA_DATA = pickle.load(open(IFILE['TCGA'], 'rb'))

	BCXX: pd.DataFrame
	TCGA: pd.DataFrame

	BCXX = BCXX_DATA['X']
	TCGA = TCGA_DATA['X']

	# Restrict to those genes
	if genes :
		BCXX = BCXX.loc[genes & set(BCXX.index), :]
		TCGA = TCGA.loc[genes & set(TCGA.index), :]

	if PARAM['zscore'] :
		BCXX = BCXX.apply(zscore, axis=1)
		TCGA = TCGA.apply(zscore, axis=1)

	BCXX = BCXX.dropna()
	TCGA = TCGA.dropna()

	assert(not BCXX.empty)
	assert(not TCGA.empty)

	# Keep only the genes that are in both tables
	# Note: the tables are indexed by the gene symbol
	(TCGA, BCXX) = TCGA.align(BCXX, join='inner', axis=0)

	# TCGA has size (gene subset) x (TCGA bulk sample)
	# BCXX has size (gene subset) x (BCXX single cell sample)
	# print(TCGA.shape, BCXX.shape)

	# TCGA and BCXX may have all-zero columns on this subset of genes;
	# The following similarity matrix will only contain
	# the nonzero columns of TCGA (as rows) and
	# the nonzero columns of BCXX (as columns);
	# Record the original columns as 'I0' and 'J0'

	# Similarity matrix for (TCGA bulk sample) x (BCXX single cell sample)
	M = cosine_similarity(TCGA, BCXX)

	return { 'M' : M, 'meta' : meta, 'I0': TCGA.columns, 'J0': BCXX.columns }


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

