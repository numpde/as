
# RA, 2017-12-04

### IMPORTS -- #

import re
import scipy
import pickle
import numpy as np
import matplotlib.pyplot as plt
from progressbar import ProgressBar as Progress

from scipy import stats

### INPUT ---- #

input_file_selected = "OUTPUT/0_selected/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

### OUTPUT --- #

extensions = ['png']
output_file_fig = "OUTPUT/2_stats_b/nnz-median.{extension}"

### MEAT ----- #

# Load the data
data = pickle.load(open(input_file_selected, "rb"))
#print(data.keys())

# Gene expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# Number of samples, number of genes
#(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# BCXX[LN][_Re]_XX
header = data['header']
#print(header)

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
	S = [(s, h) for (s, h) in enumerate(header) if re.match(g + "_[0-9]+", h)]

	# Extract the those samples from the data matrix
	GROUP[g] = {
		'data' : np.take(X, [s for (s, h) in S], axis = axis_smpl),
		'head' : [h for (s, h) in S],
		'clas' : [g] * len(S)
	}

#

print("Computing...")

ZMH = dict()
for (g, DH) in Progress()(GROUP.items()) :
	Y = DH['data']
	
	ZMH[g] = { 'nnz' : [], 'med' : [], 'hed' : DH['head'] }
	for s in range(Y.shape[axis_smpl]) :
		x = np.take(Y, s, axis=axis_smpl)
		
		ZMH[g]['nnz'].append(np.log10(np.count_nonzero(x)))
		ZMH[g]['med'].append(np.median(np.log10(x[x != 0])))

print("Plotting...")

# Each group has its own color
colors = plt.get_cmap('hsv')(np.linspace(0, 1.0, len(ZMH))).tolist()
#
for g in sorted(ZMH.keys()) :
	zmh = ZMH[g]
	
	(x, y) = (zmh['nnz'], zmh['med'])
	plt.plot(x, y, '.', color=colors.pop())
	
	# Annotate outliers
	zmh = zip(zmh['nnz'], zmh['med'], zmh['hed'])
	zmh = [(x, y, h) for (x, y, h) in zmh if (y < 0) or (y > 1.5)]
	for (x, y, h) in zmh : plt.text(x, y, h, fontsize=5)

plt.xlabel("log10(# of expressed genes)")
plt.ylabel("Median log10(expression)")
plt.legend(sorted(ZMH.keys()), loc='upper left', prop={'size': 6})

print("Saving figure...")

for ext in extensions :
	plt.savefig(output_file_fig.format(extension=ext))

