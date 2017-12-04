
# RA, 2017-11-28

### IMPORTS -- #

import pickle
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

### INPUT ---- #

input_file_selected = "OUTPUT/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

### OUTPUT --- #

output_file_plot = "OUTPUT/1_plots/silhouette_z={}.png" #.format(zscore)

### MEAT ----- #

# Load the data
data = pickle.load(open(input_file_selected, "rb"))
#
# Expression matrix
X = data['X']


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
	Y = (stats.mstats.zscore(X, axis=data['axis_smpl']) if zscore else X);
	
	# Cluster membership of each sample (BCXX)
	sample2cluster = [n[0:4] for n in data['header']]
	
	# Do the silhouette plot
	plt.clf()
	plot_silhouette(Y, sample2cluster)
	
	# Save figure
	filename = output_file_plot.format(zscore)
	plt.savefig(filename)
	print("Wrote figure to", filename)
