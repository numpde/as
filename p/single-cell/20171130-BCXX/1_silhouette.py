
# RA, 2017-11-28

### IMPORTS -- #

import re
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

### INPUT ---- #

input_file_raw_BC = "ORIGINALS/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"

### OUTPUT --- #

output_file_plot = "OUTPUT/1_silhouette_z={}.png" #.format(zscore)

### MEAT ----- #

#X = np.genfromtxt("test.txt", names=True, delimiter='\t', dtype=None)
X = np.genfromtxt(input_file_raw_BC, names=True, delimiter='\t', dtype=None)

# Column names from the raw file header
col_names = list(X.dtype.names)
#print(col_names[0:20])

# The gene type of each gene
gene_type = [t.decode("utf-8") for t in X['gene_type']]

# Histogram of gene types
gene_type_freq = set((t, len([g for g in gene_type if (g == t)])) for t in set(gene_type))
print("Gene types:", set(gene_type_freq))
# {('misc_RNA', 2034), ('TR_V_gene', 97), ('rRNA', 527), ('IG_V_pseudogene', 187), ('Mt_tRNA', 22), ('sense_overlapping', 202), ('IG_J_pseudogene', 3), ('pseudogene', 13931), ('SPIKE_IN', 3), ('sense_intronic', 742), ('TR_C_gene', 5), ('IG_D_gene', 37), ('TR_V_pseudogene', 27), ('miRNA', 3055), ('TR_J_pseudogene', 4), ('TR_D_gene', 3), ('snRNA', 1916), ('polymorphic_pseudogene', 45), ('IG_V_gene', 138), ('TR_J_gene', 74), ('lincRNA', 7114), ('antisense', 5276), ('IG_C_pseudogene', 9), ('snoRNA', 1457), ('3prime_overlapping_ncrna', 21), ('IG_J_gene', 18), ('protein_coding', 20345), ('IG_C_gene', 14), ('processed_transcript', 515), ('ERCC', 92), ('Mt_rRNA', 2)}

# Filter the BCXX(LN)_XX columns
p = re.compile(r"BC[0-9]+.*_[0-9]+", re.IGNORECASE)
col_names = [n for n in col_names if p.match(n)]

# Make a numpy matrix: genes x (metadata + samples)
X = np.asmatrix(list(X[n] for n in col_names))

# Omit genes that are not "real genes"
#
# Which genes to keep
T = {"protein_coding"}
omit = [j for (j, t) in enumerate(gene_type) if not (t in T)]
#
X = np.delete(X, omit, 1)
#
print("Omitted {} non protein coding genes".format(len(omit)))

# Omit genes that are not expressed
#
# Standard deviation over samples for each gene
std = np.std(X, axis=0).reshape(-1, 1)
# Which genes to omit?
omit = [j for j in range(0, X.shape[1]) if (std[j] == 0)]
#
X = np.delete(X, omit, 1)
#
print("Omitted {} non expressed genes".format(len(omit)))

# Size of the remaining dataset
nsamples = X.shape[0]
ngenes   = X.shape[1]
#
print("Got data with {} samples and {} nontrivial genes".format(nsamples, ngenes))

# Sanity-check
assert(len(col_names) == nsamples)
assert(nsamples), "Didn't get any samples"
assert(ngenes  ), "Didn't get any features"

# Silhouette-plot routine
def silhouette(X, C) :

	# PLOTTING CODE ADAPTED FROM
	# http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
	
	from sklearn.metrics import silhouette_samples
	
	sample_silh = silhouette_samples(X, C, metric="euclidean")
	#print(len(sample_silh))

	import matplotlib.pyplot as plt
	import matplotlib.cm as cm

	(fig, ax1) = plt.subplots(1, 1)
	fig.set_size_inches(18, 7)

	n_clusters = len(set(C))

	y_lower = 10
	for (i, c) in enumerate(reversed(sorted(list(set(C))))) :
		# Aggregate the silhouette scores for samples belonging to cluster i, and sort them
		cluster_silh = np.asarray([sv for (sv, cc) in zip(sample_silh, C) if (cc == c)])
		cluster_silh.sort()
		
		#print(i, c, len(cluster_silh), cluster_silh)
		
		size_cluster_i = cluster_silh.shape[0]
		y_upper = y_lower + size_cluster_i

		color = cm.spectral(float(i) / n_clusters)
		ax1.fill_betweenx(np.arange(y_lower, y_upper),
							0, cluster_silh,
							facecolor=color, edgecolor=color, alpha=0.7)

		# Label the silhouette plots with their cluster numbers at the middle
		ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, c)

		# Compute the new y_lower for next plot
		y_lower = y_upper + 10  # 10 for the 0 samples

# Cluster membership of each sample (BCXX)
sample2cluster = [n[0:4] for n in col_names]

# With and without the z-score transform gene-wise
for zscore in [True, False] :
	# z-score- / un- transformed data
	Y = (stats.mstats.zscore(X, axis=0) if zscore else X);
	
	# Do the silhouette plot
	plt.clf()
	silhouette(Y, sample2cluster)
	
	# Save figure
	filename = output_file_plot.format(zscore)
	plt.savefig(filename)
	print("Wrote figure to", filename)
