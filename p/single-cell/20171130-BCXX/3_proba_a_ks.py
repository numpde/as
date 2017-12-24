
# RA, 2017-12-05 (Initial)
# RA, 2017-12-21 (Revision)

## ================== IMPORTS :

import os, re, sys, math, pickle, inspect

import numpy as np
import matplotlib.pyplot as plt

from scipy           import stats
from multiprocessing import cpu_count
from joblib          import Parallel, delayed
from progressbar     import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	'BC data' : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), ("File not found: " + f)

## =================== OUTPUT :

OFILE = {
	'gene-ks' : "OUTPUT/3_proba_a/UV/gene-ks.pkl",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## =================== PARAMS :

PARAM = {
	# Number of parallel computing processes
	'#jobs' : min(10, math.ceil(3/4 * cpu_count())),
	
	# Clusters/groups by 'p'atient or by 'b'atch?
	'groups' : 'b',
	
	# Suspicious samples to be omitted
	'omit samples' : {
		"BC07LN_20", 
	},
}


# Test mode
TESTMODE = ("TEST" in sys.argv)

## ==================== PREPA :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

## ====================== (!) :

# (alphabetical order)
KS_meta = ['max', 'mean', 'median', 'min', 'std']

# Differential expression within a collection of empirical proba
def KS(P) :
	ks = [
		stats.ks_2samp(p, q)[0] 
		for (m, p) in enumerate(P) 
		for (n, q) in enumerate(P) if (m < n)
	]
	# (alphabetical order)
	return (np.max(ks), np.mean(ks), np.median(ks), np.min(ks), np.std(ks))

## ===================== WORK :

# Load the BC data
data = pickle.load(open(IFILE['BC data'], "rb"))
print("Data keys:", list(data.keys()))

# Expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# Labels of samples of the form BCXX[LN][_Re]_XX 
sample_labels = data['header']

# Remove strange samples
for h in PARAM['omit samples'] :
	X = np.delete(X, sample_labels.index(h), axis=axis_smpl)
	sample_labels.remove(h)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# Get the sample groups
#
if (PARAM['groups'] == 'p') : 
	# Option 1: BCXX -- by patient
	groups = [n[0:4] for n in sample_labels]
#
if (PARAM['groups'] == 'b') : 
	# Option 2: BCXX[LN][_Re] -- by batch
	groups = [re.findall("(.*)_[0-9]+", n)[0] for n in sample_labels]
#
assert(groups), "Grouping failed"
#
# Make group labels unique and sort
groups = sorted(list(set(groups)))
print("Groups:", ', '.join(groups))

# Collect sample labels of the form "g"_XX by group
G2H = {
	g : tuple(h for h in sample_labels if re.match(g + "_[0-9]+", h))
	for g in groups 
}

# Collect the IDs of those samples
G2S = {
	g : tuple(sample_labels.index(h) for h in H)
	for (g, H) in G2H.items()
}

# Split the expression matrix by groups, arrange by genes
G2X = {
	g : np.moveaxis(np.take(X, S, axis=axis_smpl), axis_gene, 0)
	for (g, S) in G2S.items()
}

# Compute differential expression for gene #n
def job(n) :
	return (data['gene_id'][n], KS([G2X[g][n] for g in groups]))

# E2KS : ENSG --> (Gene differential expression)
E2KS = dict(
	Parallel(n_jobs=PARAM['#jobs'])(
		delayed(job)(n) 
		for n in Progress()(range(n_genes))
	)
)

if TESTMODE : exit()

# Save results
pickle.dump(
	{ 
		'E2KS'    : E2KS, 
		'KS_meta' : KS_meta,
		'groups'  : groups,
		'G2H'     : G2H,
		'G2S'     : G2S,
		'PARAM'   : PARAM,
		'script'  : THIS,
	}, 
	open(OFILE['gene-ks'], "wb")
)
