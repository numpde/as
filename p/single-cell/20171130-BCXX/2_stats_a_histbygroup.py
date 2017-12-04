
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

output_file_fig = "OUTPUT/2_stats_a/{normalization}_{group}_XX.{extension}"

### MEAT ----- #

# Load the data
data = pickle.load(open(input_file_selected, "rb"))
#print(data.keys())

# Expression matrix
X = data['X']

# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# Number of samples, number of genes
#(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])

# BCXX(LN)(_Re)_XX
header = data['header']
#print(header)

# Get the sample groups
#
# Option 1: BCXX
#groups = [n[0:4] for n in header]
#
# Option 2: BCXX(LN)(_Re)
groups = [re.findall("(.*)_[0-9]+", n)[0] for n in header]

# Make groups unique and sort
groups = sorted(list(set(groups)))
print("Groups:", ', '.join(groups))

GROUP = dict()
for g in groups :
	# Will omit BCXXLN_XX and BCXX_Re_XX samples
	# by matching only BCXX_XX for BCXX = g
	S = [s for (s, h) in enumerate(header) if re.match(g + "_[0-9]+", h)]

	# Extract the those samples from the data matrix
	GROUP[g] = np.take(X, S, axis = axis_smpl)


# https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
from scipy.stats import gaussian_kde

for normalize in [True, False] :
	print("Computing with normalization flag:", normalize)

	xlabel = ("log10(normalized expression)" if normalize else "log10(expression)")
	
	TF = dict()
	for (g, Y) in Progress()(GROUP.items()) :
		TF[g] = []
		for s in range(Y.shape[axis_smpl]) :
			x = np.take(Y, s, axis = axis_smpl)
		
			x = x[x != 0]

			d = dict()
			d['median'] = np.median(x)
			
			if normalize : x = x / np.sum(x)
			
			x = np.log10(x)
			
			t = np.linspace(min(x), max(x), 100)
			f = gaussian_kde(x)(t)
			
			TF[g].append( (t, f, d) )

	print("Plotting...")

	mint = min(min(t) for tf in TF.values() for (t, _, _) in tf)
	maxt = max(max(t) for tf in TF.values() for (t, _, _) in tf)
	maxf = max(max(f) for tf in TF.values() for (_, f, _) in tf)

	for (g, tfd) in Progress()(TF.items()) :
		plt.clf()
		
		for (t, f, d) in tfd :
			is_alive = (d['median'] > 1)
			plt.plot(t, f, ('C2' if is_alive else 'C1'))
		
		plt.xlim(( mint, maxt ))
		plt.ylim((    0, maxf ))
			
		plt.title("Sample subset: {}_XX ({} samples)".format(g, Y.shape[axis_smpl]))
		plt.xlabel(xlabel)
		plt.ylabel("Frequency")

		plt.savefig(
			output_file_fig.format(
				normalization=("sum1" if normalize else "orig"), 
				group=g, 
				extension="png"
			)
		)
