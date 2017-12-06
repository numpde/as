
# RA, 2017-12-05

### IMPORTS -- #

import re
import math
import pickle
import inspect
import numpy as np
import matplotlib.pyplot as plt

from scipy           import stats
from multiprocessing import cpu_count
from joblib          import Parallel, delayed
from progressbar     import ProgressBar as Progress

### INPUT ---- #

input_file_selected = "OUTPUT/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

### OUTPUT --- #

output_file_e2ks = "OUTPUT/3_proba_a/e2ks.pkl"

### PARAM ---- #

# Number of parallel computing processes
num_procs = math.ceil(cpu_count() / 2)

### (!!!) ---- #


# Differential expression between two collections of empirical proba
def DE(P, Q) :
	# Distance between two empirical proba
	def dist(p, q) : return stats.ks_2samp(p, q)[0]

	return np.max([dist(p, q) for p in P for q in Q])

### MEAT ----- #

# Load the BC datas
data = pickle.load(open(input_file_selected, "rb"))
#print(data.keys())

# Expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# Labels of samples of the form BCXX[LN][_Re]_XX 
header = data['header']
#print(header)

# Remove strange samples
for h in { "BC07LN_20" } :
	X = np.delete(X, header.index(h), axis=axis_smpl)
	header.remove(h)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# Get the sample groups
#
# Option 1: BCXX
groups = [n[0:4] for n in header]
#
# Option 2: BCXX[LN][_Re]
#groups = [re.findall("(.*)_[0-9]+", n)[0] for n in header]

# Make groups unique and sort
groups = sorted(list(set(groups)))
print("Groups:", ', '.join(groups))

# Collect numbers of samples of the form "g"_XX by group
G2S = {
	g : [s for (s, h) in enumerate(header) if re.match(g + "_[0-9]+", h)]
	for g in groups 
}

# Split the expression matrix by groups
G2X = {
	g : np.take(X, G2S[g], axis=axis_smpl) 
	for g in groups 
}

# Compute differential expression for gene #n
def job(n) :
	P = [np.take(G2X[g], n, axis=axis_gene) for g in groups]
	ks = DE(P, P)
	
	return (data['gene_id'][n], ks)

# [Gene Ensembl ID] --> [Gene differential expression]
E2KS = dict(Parallel(n_jobs=num_procs)(delayed(job)(n) for n in Progress()(range(n_genes))))

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
script = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Save results
pickle.dump({ 'E2DE' : E2KS, 'script' : script }, open(output_file_e2ks, "wb"))
