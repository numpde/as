
# RA, 2017-11-28

### IMPORTS -- #

import re
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

### INPUT ---- #

input_file_raw_BC = "ORIGINALS/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
input_file_select = "ORIGINALS/gene-selection-txp.txt"

### OUTPUT --- #

output_file_plot = "OUTPUT/1_silhouette_z={}.png" #.format(zscore)

### MEAT ----- #

# Load the raw data
X = np.genfromtxt(input_file_raw_BC, names=True, delimiter='\t', dtype=None)

# Feature tags from the raw file 
header = list(X.dtype.names)
#print(header[0:20])

# Filter the BCXX(LN)_XX columns
p = re.compile(r"BC[0-9]+.*_[0-9]+", re.IGNORECASE)
header = [n for n in header if p.match(n)]

# The gene type of each gene
gene_type = [t.decode("utf-8") for t in X['gene_type']]

# The ENSGXXXXXXXXXXX gene id of each gene
gene_id = [t.decode("utf-8")[0:15] for t in X['gene_id']]

# Histogram of gene types
#gene_type_freq = set((t, len([g for g in gene_type if (g == t)])) for t in set(gene_type))
#print("Gene types:", set(gene_type_freq))
## {('misc_RNA', 2034), ('TR_V_gene', 97), ('rRNA', 527), ('IG_V_pseudogene', 187), ('Mt_tRNA', 22), ('sense_overlapping', 202), ('IG_J_pseudogene', 3), ('pseudogene', 13931), ('SPIKE_IN', 3), ('sense_intronic', 742), ('TR_C_gene', 5), ('IG_D_gene', 37), ('TR_V_pseudogene', 27), ('miRNA', 3055), ('TR_J_pseudogene', 4), ('TR_D_gene', 3), ('snRNA', 1916), ('polymorphic_pseudogene', 45), ('IG_V_gene', 138), ('TR_J_gene', 74), ('lincRNA', 7114), ('antisense', 5276), ('IG_C_pseudogene', 9), ('snoRNA', 1457), ('3prime_overlapping_ncrna', 21), ('IG_J_gene', 18), ('protein_coding', 20345), ('IG_C_gene', 14), ('processed_transcript', 515), ('ERCC', 92), ('Mt_rRNA', 2)}

# Make a numpy matrix: (metadata + samples) x genes
X = np.asarray(list(X[n] for n in header))

# Along which axis/dimension of the data are the samples / the genes
(axis_sample, axis_gene) = (0, 1)

# Cluster membership of each sample (BCXX)
sample2cluster = [n[0:4] for n in header]

if False :
	# Omit the BC01 cluster
	I = [i for (i, c) in enumerate(sample2cluster) if (c == "BC01")]
	X = np.delete(X, I, axis_sample)
	sample2cluster = [c for (i, c) in enumerate(sample2cluster) if (i not in I)]
	del header # invalidate inconsitent variable

# Omit irrelevant genes 

if False : 
	# Keep only genes of this type:
	T = {"protein_coding"}
	
	# Which genes to omit?
	J = [j for (j, t) in enumerate(gene_type) if (t not in T)]
	X = np.delete(X, J, axis_gene)
	print("Omitted {} genes that are not {}".format(len(J), T))

if True :
	# Use the selection file to filter relevant genes
	Y = np.genfromtxt(input_file_select, names=True, delimiter='\t', dtype=None)
	
	# 15 is the length of the ENSGXXXXXXXXXXX gene_id code
	Y = set([t.decode("utf-8")[0:15] for t in Y['gene_id']])
	
	# Which genes to omit?
	J = [j for (j, i) in enumerate(gene_id) if not (i in Y)]
	X = np.delete(X, J, axis_gene)
	print("Omitted {} genes not listed in {}".format(len(J), input_file_select))

if True :
	# Omit genes that are not expressed
	J = (np.std(X, axis=axis_sample) == 0)
	X = X[:, np.invert(J)]
	print("Omitted {} non-expressed genes".format(sum(J)))

# Size of the remaining dataset
(n_samples, n_genes) = (X.shape[axis_sample], X.shape[axis_gene])
#
print("Got data with {} samples and {} nontrivial genes".format(n_samples, n_genes))

# Sanity check
assert(n_samples), "Didn't get any samples"
assert(n_genes  ), "Didn't get any genes"

# Silhouette-plot routine
# X = data with samples in rows
# C = sample cluster labels
def plot_silhouette(X, C) :

	# PLOTTING CODE ADAPTED FROM
	# http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
	
	assert(len(C) == X.shape[0]), "C-X dimensions inconsistent"
	
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
		ax1.fill_betweenx(
			np.arange(y_lower, y_upper),
			0, cluster_silh,
			facecolor=color, edgecolor=color, alpha=0.7)

		# Label the silhouette plots with their cluster numbers at the middle
		ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, c)

		# Compute the new y_lower for next plot
		y_lower = y_upper + 10  # 10 for the 0 samples


# With and without gene-wise z-score transform
for zscore in [True, False] :
	
	# Get the z-score-transformed or the original data
	# If transform then gene-wise, i.e. over the sample axis/dimension
	Y = (stats.mstats.zscore(X, axis=axis_sample) if zscore else X);
	
	# Do the silhouette plot
	plt.clf()
	plot_silhouette(Y, sample2cluster)
	
	# Save figure
	filename = output_file_plot.format(zscore)
	plt.savefig(filename)
	print("Wrote figure to", filename)
