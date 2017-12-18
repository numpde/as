
# RA, 2017-11-28

### IMPORTS -- #

import re
import pickle
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

### INPUT ---- #

input_file_selected = "OUTPUT/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

### OUTPUT --- #

pass

### MEAT ----- #

# Load the data
data = pickle.load(open(input_file_selected, "rb"))
#
# Expression matrix
X = data['X']


# Labels for axis/dimension of data
(axis_smpl, axis_gene) = (data['axis_smpl'], data['axis_gene'])

# Labels of samples of the form BCXX[LN][_Re]_XX 
sample_labels = data['header']

# Remove strange samples
for h in { "BC07LN_20" } :
	X = np.delete(X, sample_labels.index(h), axis=axis_smpl)
	sample_labels.remove(h)

# Number of samples / genes in the expression matrix
(n_samples, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])



# Get the sample groups
#
# Option 1: BCXX -- by patient
#groups = [n[0:4] for n in sample_labels]
#
# Option 2: BCXX[LN][_Re] -- by batch
groups = [re.findall("(.*)_[0-9]+", n)[0] for n in sample_labels]

# Make groups unique and sort
groups = sorted(list(set(groups)))
print("Groups:", ', '.join(groups))

# Collect numbers of samples of the form "g"_XX by group
G2S = {
	g : [s for (s, h) in enumerate(sample_labels) if re.match(g + "_[0-9]+", h)]
	for g in groups 
}


#print(groups)

#plt.ion()

gene3 = np.take(X, 3, axis=axis_gene)
for (n, g) in enumerate(sorted(groups)) :
	bc = np.take(gene3, G2S[g], axis=axis_smpl)
	#print(g, bc)
	
	plt.hist(bc, histtype='step', log=True)
	
plt.show()

#input("")
