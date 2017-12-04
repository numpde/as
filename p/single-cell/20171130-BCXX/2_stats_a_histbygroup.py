
# RA, 2017-12-02

### IMPORTS -- #

import re
import scipy
import pickle
import numpy as np
import matplotlib.pyplot as plt
from progressbar import ProgressBar as Progress

from scipy import stats

### INPUT ---- #

input_file_selected = "OUTPUT/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

### OUTPUT --- #

output_file_fig = "OUTPUT/2_stats_b/{}.png" #.format(group)

### MEAT ----- #

# Load the data
data = pickle.load(open(input_file_selected, "rb"))
#print(data.keys())

# Expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

#
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# BCXX(LN)(_Re)_XX
header = data['header']
#print(header)

# BCXX
#groups = [n[0:4] for n in header]

# BCXX(LN)(_Re)
groups = [re.findall("(.*)_[0-9]+", n)[0] for n in header]

# Make unique and sort
groups = list(sorted(set(groups)))
print(groups)

GROUP = dict()
for g in groups :
	# Will omit BCXXLN_XX and BCXX_Re_XX samples
	# by matching only BCXX_XX for BCXX = g
	S = [s for (s, h) in enumerate(header) if re.match(g + "_[0-9]+", h)]

	# Extract the those samples from the data matrix
	GROUP[g] = np.take(X, S, axis = axis_smpl)


# https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
from scipy.stats import gaussian_kde

for (g, Y) in Progress()(GROUP.items()) :
	plt.clf()
	plt.title("Sample subset: {}_XX ({} samples)".format(g, Y.shape[axis_smpl]))
	plt.xlabel("log10(expression)")
	plt.ylabel("Frequency")
	
	for s in range(Y.shape[axis_smpl]) :
		x = np.take(Y, s, axis = axis_smpl)
		
		#print(len(x), np.histogram(np.log(x[x != 0])))
		
		x = np.log10(x[x != 0])
		
		t = np.linspace(min(x), max(x), 100)
		f = gaussian_kde(x)(t)
		
		plt.plot(t, f)

	plt.savefig(output_file_fig.format(g))

