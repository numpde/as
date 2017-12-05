
# RA, 2017-12-02

# (In progress)

### IMPORTS -- #

import re
import scipy
import pickle
import inspect
import numpy as np
import matplotlib.pyplot as plt
from progressbar import ProgressBar as Progress

from scipy import stats

### INPUT ---- #

input_file_selected = "OUTPUT/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

### OUTPUT --- #

output_file_e2ks = "OUTPUT/3_proba_a/e2ks.pkl"

### MEAT ----- #

# Load the BC data
data = pickle.load(open(input_file_selected, "rb"))
#print(data.keys())

# Expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# BCXX(LN)(_Re)_XX
header = data['header']
#print(header)

# Remove strange samples
strange = { "BC07LN_20" }
for h in strange :
	s = header.index(h)
	header.pop(s)
	X = np.delete(X, s, axis=axis_smpl)


#
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# Get the sample groups
#
# Option 1: BCXX
#groups = [n[0:4] for n in header]
#
# Option 2: BCXX[LN][_Re]
groups = [re.findall("(.*)_[0-9]+", n)[0] for n in header]

# Make groups unique and sort
groups = sorted(list(set(groups)))
print("Groups:", ', '.join(groups))

GROUP = dict()
for g in groups :
	# Find samples of the form [g]_XX
	S = [s for (s, h) in enumerate(header) if re.match(g + "_[0-9]+", h)]

	# Extract the those samples from the data matrix
	GROUP[g] = np.take(X, S, axis = axis_smpl)


E2KS = dict()

for n in Progress()(range(n_genes)) :
	
	ks = np.zeros( (len(groups), len(groups)) )
	def f(Y) : return np.take(Y, n, axis = axis_gene)
	for (i, a) in enumerate(groups) :
		for (j, b) in enumerate(groups) :
			A = f(GROUP[a]); A = A[A.nonzero()]
			B = f(GROUP[b]); B = B[B.nonzero()]
			if len(A) and len(B) :
				(d, p) = scipy.stats.ks_2samp(A, B)
			else :
				d = 0
			ks[i, j] += d
	
	E2KS[data['gene_id'][n]] = np.sum(ks)
	
	#if (n > 100) :
		#print("Warning: early abort.")
		#break


# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
script = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Save results
pickle.dump({ 'E2KS' : E2KS, 'script' : script }, open(output_file_e2ks, "wb"))
