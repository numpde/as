
# RA, 2018-01-30

# Combines g, h1 and h2 ?
# Outdated compared to those.


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
	# SIMILARITY COMPUTATION
	
	# o) Expression data
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
	
	# o) Gene sets
	'genesets' : "OUTPUT/f_gene_subsets/subsets.pkl",
	
	# CLASSIFICATION
	
	# o) TCGA PAM50 data
	'PAM50' : "OUTPUT/e_prepared/UV/tcga.pkl",
	
	# o) TCGA clinical data
	'clinical' : "OUTPUT/b1_clinical/tcga-clinical.pkl",
}


## =================== OUTPUT :

OFILE = {
	'classified' : "OUTPUT/h3_class/{ext}/{A}-{B}/{method}/{geneset}.pdf",
}


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.2)),
	
	# Record just in case
	'testmode' : TESTMODE,
	
	# Do columns-wise normalization?
	'normcol' : False,
	
	# Do gene-wise zscore transform?
	'zscore' : True,
	
	# Classification methods
	'Classification' : {
		'PAM50' : {
			'key' : 'PAM50',
			'sub' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
		},
		
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
	},
	
	# Extension for plots
	'ext' : { 'pdf', 'png' },
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# zscore transform of a pandas series
def zscore(S) : return (S - np.mean(S)) / np.std(S)

# normalization of a pandas series
def normcol(S) : return S / np.sum(np.abs(S))

# softmax of a pandas series
def softmax(S, amp=1) : return np.exp(amp*S) / np.sum(np.exp(amp*S))

# Change GO:000 to GO-000 in filenames
def nicer(filename) : return filename.replace("GO:", "GO-")

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
		# Guess that L is an iterable
		return [abbr(i) for i in L]


## ====================== (!) :

# 
def compile_class_table() :
	
	class_keys = [p['key'] for (m, p) in PARAM['Classification'].items()]
	
	# Clinical
	C = pickle.load(open(IFILE['clinical'], 'rb'))['patient_barcode']
	# Identify record by patient barcode
	assert(C.columns.is_unique)
	
	# PAM50
	P = pickle.load(open(IFILE['PAM50'], 'rb'))['subtype']
	# Identify record by aliquot barcode
	assert(P['aliquot_barcode'].is_unique)
	# Some patients have multiple aliquots though
	assert(not P['patient_barcode'].is_unique)
	# Use the aliquot barcode as index
	P = P.set_index('aliquot_barcode')
	
	for p in P['patient_barcode'] :
		tt = P['tissue_type'][P['patient_barcode'] == p]
		assert((tt == 'tumor').sum() == 1)
	
	# TCGA aliquots
	A = pickle.load(open(IFILE['TCGA'], 'rb'))['X'].columns
	
	# Merged classification table
	M = pd.DataFrame(index=class_keys, columns=A)
	
	# Collect available classification for each aliquot
	for a in M.columns :
		if (a in P.index) :
			M.loc[ 'PAM50', a ] = P.loc[ a, 'PAM50' ]
		
		# Skip aliquot if we know it is not from tumor ...
		
		if (a in P.index) :
			if (P['tissue_type'][a] != 'tumor') : 
				continue
		
		# ... otherwise extract clinical info by patient barcode
		
		# Patient barcode
		p = a[0:12]
		
		if (p in C.columns) :
			for k in M.index :
				if (k in C.index) :
					M[a][k] = C[p][k]
	
	return M


## ====================== (!) :


# Cosine similarity between two pandas dataframes (column-wise)
# https://en.wikipedia.org/wiki/Cosine_similarity
def cosine_similarity(A, B) :
	
	# Drop zero columns
	A = A.loc[ :, (A != 0).any(axis=0) ]
	B = B.loc[ :, (B != 0).any(axis=0) ]
	
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


# Compute the similarity along a subset of genes
def compute_similarity(A, B, genes=None) : 
	
	# Local copy
	A = A.copy()
	B = B.copy()
	
	# Restrict to those genes
	if genes :
		A = A.loc[genes]
		B = B.loc[genes]
	
	# Apply columns-wise normalization
	if PARAM['normcol'] :
		A = A.apply(normcol, axis=0)
		B = B.apply(normcol, axis=0)
	
	A = A.dropna()
	B = B.dropna()
	
	# Apply zscore transform gene-wise
	if PARAM['zscore'] :
		A = A.apply(zscore, axis=1)
		B = B.apply(zscore, axis=1)
	
	A = A.dropna()
	B = B.dropna()
	
	# Keep only the genes that are in both tables
	# Note: the tables are indexed by the gene symbol
	(A, B) = A.align(B, join='inner', axis=0)
	
	# Similarity matrix for (TCGA bulk sample) x (BCXX single cell sample)
	M = cosine_similarity(A, B)
	
	return M


## ===================== WORK :

# Compute the similarity grouping by class
def classify(b, genes, method) :
	
	assert(b in ['TCGA', 'BCXX'])
	
	# TCGA subtype classification data
	C = compile_class_table()
	
	# Load TCGA expression data
	A = pickle.load(open(IFILE['TCGA'], 'rb'))['X']
	
	# Keep only the TCGA samples that have been classified
	A = A.loc[ :, list(pd.Series(A.columns).isin(C.columns)) ]
	
	if (b == 'TCGA') : B = A.copy()
	if (b == 'BCXX') : B = pickle.load(open(IFILE['BCXX'], 'rb'))['X']
	
	# 
	M = compute_similarity(A, B, genes)
	
	# Group into a (TCGA subtype) x (BCXX tumor) matrix
	def group(M) :
		
		# Omit records with missing data
		N = M.dropna().astype(float)
		
		# Group by class
		N = N.groupby(
			# Aliquot barcode to class
			lambda a : C[a][method['key']],
			axis=0
		).mean()
		
		# Use the order of the subtypes list
		N = N.reindex(pd.Categorical(N.index, categories=method['sub']).sort_values(), axis=0)
		
		assert(b in ['TCGA', 'BCXX'])
	
		if (b == 'BCXX') :
			# Group by BCXX tumors
			N = N.groupby(
				# Sample "BCXX[LN][_Re]_YY" to cluster "BCXX"
				lambda c : c[:-3],
				axis=1
			).mean()
		
		if (b == 'TCGA') :
			# Group by class
			N = N.groupby(
				# Aliquot barcode to class
				lambda a : C[a][method['key']],
				axis=1
			).mean()
			
			N = N.reindex(pd.Categorical(N.columns, categories=method['sub']).sort_values(ascending=False), axis=1)
		
		return N
	
	# Group into a (TCGA subtype) x (BCXX tumor) matrix
	N = group(M)

	return N

# Plot a (TCGA subtype) x (BCXX tumor) matrix
def plot(N, m) :
	
	# Entries are potentially in [-1, 1]
	for c in N.columns : N[c] = softmax(N[c], amp=3)
	
	## Normalize row-wise
	#for i in N.index : N.loc[i] /= N.loc[i].sum()
	
	plt.clf()
	im = plt.imshow(N.as_matrix(), aspect='auto', cmap=plt.cm.Blues, vmin=0, vmax=1)
	plt.xticks(range(len(N.columns)), abbr(list(N.columns)), rotation=60)
	plt.yticks(range(len(N.index)), abbr(list(N.index)), rotation=60)
	if (len(N.columns) > 6) :
		plt.tick_params(axis='x', which='major', labelsize=5)
	#plt.tick_params(axis='y', which='major', labelsize=10)
	#plt.colorbar()
	
	plt.ylabel(m + " normed cos-sim to...")
	ax = plt.gca()
	
	# https://matplotlib.org/users/tight_layout_guide.html
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	cax = make_axes_locatable(plt.gca()).append_axes("right", "5%", pad="3%")
	plt.colorbar(im, cax=cax)
	plt.tight_layout()
	cax.tick_params(labelsize=5) 

	
	plt.sca(ax)

# 
def compute(b, geneset) :
	
	for (m, method) in PARAM['Classification'].items() :
		
		#print("Gene set:", setid)
		#print("Method:", method)
	
		N = classify(b, geneset['set'], method)
		
		plot(N, m)
		plt.title("{} ({})".format(geneset['info'], len(geneset['set'])), fontsize=6)
		#plt.show()
		
		for ext in PARAM['ext'] : 
			# Filename for figure
			f = OFILE['classified'].format(ext=ext, A='TCGA', B=b, method=m, geneset=nicer(geneset['id']))
			# Create output directories
			os.makedirs(os.path.dirname(f), exist_ok=True)
			# Save figure
			plt.savefig(f)

def COMPUTE() :
	
	# Load gene subsets of interest
	gene_subsets = pickle.load(open(IFILE['genesets'], 'rb'))['S']
	
	for b in ['TCGA', 'BCXX'] : 
		print("Classifying", b)
		Parallel(n_jobs=PARAM['#proc'])(
			delayed(compute)(b, geneset)
			for geneset in Progress()(gene_subsets)
		)


## ==================== ENTRY :

if (__name__ == "__main__") :
	COMPUTE()

