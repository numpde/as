
# RA, 2018-01-30

# Run as
#    python3 i2*.py


## ============ LOCAL IMPORTS :

from utils import make_keras_picklable


## ================== IMPORTS :

import os
import sys
import math
import pickle
import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from itertools import chain

from keras import backend as keras_backend


## ==================== INPUT :

IFILE = {
	# Expression data
	#'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
	
	# Trained keras model
	'predictor' : "OUTPUT/i1_pam50/UV/trained.pkl",
}


## =================== OUTPUT :

OFILE = {
	'pred-sc-samples' : "OUTPUT/i2_pam50_bcxx/pred_singlecell/samples.{ext}",
	'pred-sc-bygroup' : "OUTPUT/i2_pam50_bcxx/pred_singlecell/bygroup/{group}.{ext}",
	
	'pred-bulk' : "OUTPUT/i2_pam50_bcxx/pred_bulk/bygroup/{group}.{ext}",
	
	'pred-sc2bulk-img' : "OUTPUT/i2_pam50_bcxx/pred_sc2bulk/img/conf-win={win}.{ext}",
	'pred-sc2bulk-dat' : "OUTPUT/i2_pam50_bcxx/pred_sc2bulk/dat/offendants.pkl",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	## Number of parallel computing processes
	#'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.2)),
	
	# Record just in case
	'testmode' : TESTMODE,
	
	'PAM50-types' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
	
	# Max number of cells in the artificial tumor
	'sc2bulk-N' : 20,
	
	## Extension for plots
	#'ext' : { 'pdf', 'png' },
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Check if pandas series has unique items
def is_unique(S) : return (S.unique().size == S.size)


## ====================== (!) :

# Infer tumor labels from sample labels for BCXX data
def bcxx_groups(P) :
	# Group by tumor
	G2H = dict(P.groupby(by=(lambda c : c[:-3]), axis=1).groups)
	
	# Sorted group names
	G = sorted(G2H.keys())
	
	# Sample groups in the same order
	S = [G2H[g] for g in G]
	
	return (G, S)


# Normalization as in the training
def normalize(A) :
	
	# Normalize columns-wise
	for c in A.columns : A[c] /= (sum(A[c]) or 1)
	
	# Log-transform
	A = np.log(1 + 10 * A)
	
	return A


def predict(M, B) :
	
	# Classify samples and assemble into a dataframe
	# Normalization should be done only once!
	P = pd.DataFrame(M['m'].predict(normalize(B).as_matrix().T).T, index=M['L'], columns=B.columns)

	return P
	

def predict_bcxx(singlecell=True) :
	# BCXX expression as pandas table
	B = pickle.load(open(IFILE['BCXX'], 'rb'))['X']
	
	# Trained PAM50 predictor bundle-+*
	M = pickle.load(open(IFILE['predictor'], 'rb'))
	
	# Use spike-in to reweight samples
	W = np.sqrt(B[B.index == "EC2"].values * B[B.index == "EC15"].values).flatten()
	for (w, c) in zip(W/sum(W), B.columns) : B[c] *= (w or 1)
	
	# Select genes = features
	B = B.loc[ M['F'], : ]
	
	## Omit "dying" cells
	#def drop(s) : return (s < s.median()/5)
	#B = B.loc[ :, ~drop(B.sum(axis=0)) ]
	
	# Simulate bulk sequencing?
	if not singlecell :
		# Cell ID to tumor ID
		c2t = (lambda c : c[:-3])
		# Get the average expression in each tumor
		B = B.groupby(by=c2t, axis=1).mean()
	
	# 
	P = predict(M, B)
	
	return (B, M, P)


## ===================== WORK :


## == SC-TO-BULK PREDICTIONS == ##

def plot_sc2bulk() :
	
	(B, M, P) = predict_bcxx(singlecell=True)
	
	# Class labels
	L = M['L']
	
	plt.figure(figsize=(6, 4), dpi=100)
	
	# Handles for the legend
	H = []
	
	# List of "worst offendant" artificial tumors by class
	T = dict()
	
	# Attempt to construct a tumor of subtype c0
	for c0 in L :
		# Sort the classification by their proximity to class c0
		P = P.sort_values(by=c0, axis=1, ascending=True)
		# Sort the samples accordingly
		b = B[P.columns]
		# Subtype prediction for artificial subtumors
		def assemble_tumor(c) :
			# Take the columns up to column 'c'
			t = b.loc[:, :c]
			# Take at most N cells
			return t[t.columns[-PARAM['sc2bulk-N']:]]
		#
		p = [
			predict(M, pd.DataFrame(data={c : assemble_tumor(c).mean(axis=1)}))
			for (n, c) in enumerate(b.columns)
		]
		# Convert this to a dataframe
		p = pd.concat(p, join='inner', axis=1)
		
		f = P.ix[c0].values   # single cell confidence
		n = np.arange(len(f)) # number of cells in bulk
		y = p.ix[c0].values   # bulk confidence
		
		# Cut-off below this confidence
		minf = 0.01
		(f, n, y) = (f[f >= minf], n[f >= minf], y[f >= minf])
		
		# https://matplotlib.org/users/colormaps.html
		cmap = [plt.cm.Purples, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Reds][L.index(c0)]
		
		markers = ['o', '^', 'v', '*', 's']
		
		h0 = plt.plot(f, y, linestyle='-', marker=markers[L.index(c0)], color=cmap(200), markersize=3, lw=1)[0]
		#h1 = plt.scatter(f, y, c=f, cmap=cmap, s=10)
		
		# The "worst offendant"
		m = np.argmax(y * (f <= 0.1))
		#plt.scatter(f[m], y[m], s=30, c=cmap(200))
		t = assemble_tumor(b.columns[n[m]])
		T[c0] = t
		
		H.append(h0)

	plt.legend(H, L, loc='lower left')
	
	plt.xlabel("Upper confidence bound for single cell")
	plt.ylabel("Bulk confidence for the same subtype")
	
	plt.xticks(np.linspace(0, 1, 11))
	plt.yticks(np.linspace(0, 1, 11))
	
	plt.ylim((-0.03, 1.03))
	plt.tight_layout()
	
	# Save the offendants to disk
	
	pickle.dump({'worst' : T}, open(OFILE['pred-sc2bulk-dat'], 'wb'))
	
	# Save image extracts to disk
	
	for x in [0, 1, 2, 3, 4, 5] :
		plt.xlim(((x/10) - 0.01, (0.5 + x/10) + 0.01))
		plt.gca().invert_xaxis()
		plt.savefig(OFILE['pred-sc2bulk-img'].format(win=x, ext="pdf"))
	
	plt.close()


## == SINGLE CELL PREDICTIONS == ##

# Plot the predictions of the classifier for all samples
def plot_singlecell_samples(P) :
	
	(G, S) = bcxx_groups(P)
	
	# Sample label to (predicted class, -confidence)
	def key(h) : 
		n = np.argmax(P[h].values)
		return (n, -P[h][n])
	
	# Sort sample groups by predicted class and -confidence
	S = [sorted(s, key=key) for s in S]
	
	# Widths of each cluster
	R = [len(s) for s in S]
	# Insert spacing between them
	spacing = 2
	R = list(chain(*zip([spacing] * len(R), R)))[1:]
	# Add space for the subtypes labels and the colorbar
	R = [sum(R)*0.03] + [spacing] + R + [sum(R)*0.03, sum(R)*0.01, spacing]
	R = np.cumsum(R)
	R = R / max(R)
	
	fig = plt.figure(figsize=(10, 3), dpi=300, frameon=False)
	
	# Vertical orientation
	origin = ["upper", "lower"][0]
	
	# Ploting axes for the groups
	AX = [
		fig.add_axes([a, 0, b-a, 1])
		for (a, b, _) in zip(R[1::2], R[2::2], S)
		# The third component is to have the correct number of axes
	]
	
	# Group-wise sample classification
	
	for (g, s, ax) in zip(G, S, AX) :
		
		im = ax.imshow(
			P[s].as_matrix(), 
			aspect='auto', origin=origin, cmap=plt.cm.Blues,
			vmin=0, vmax=1
		)
		
		ax.set_xticks([])
		ax.set_yticks([])
		
		# Group label
		ax.text(
			0.95, 0.99, g, 
			size=6,
			transform=ax.transAxes, rotation=90, 
			ha='right', va='top',
			color="red"
		)
	
	# Colorbar
	
	cax = fig.add_axes([R[-3], 0.05, R[-2]-R[-3], 0.9])
	cb = plt.colorbar(im, cax=cax, ticks=[0, 0.5, 1])
	#
	cb.set_clim(0, 1)
	cax.set_yticklabels(["0%", "confidence", "100%"], va='center', rotation=90, size=5)
	cax.yaxis.set_ticks_position('left')
	cax.tick_params(axis='both', which='both', length=0)
	
	# Class labels
	
	ax = fig.add_axes([0, 0, R[0], 1])
	ax.axis('off')
	#
	L = list(P.index)
	if (origin == "upper") : L = list(reversed(L))
	#
	for (y, t) in zip(np.linspace(0, 1, 1 + 2*len(L))[1::2], L) :
		
		ax.text(
			0.5, y, t, 
			size=8,
			transform=ax.transAxes, rotation=90, 
			ha='center', va='center',
			color='blue',
		)
	
	# Save to disk
	
	plt.savefig(OFILE['pred-sc-samples'].format(ext="pdf"))
	plt.close()


def plot_singlecell_subtypes(P) :
	
	(G, S) = bcxx_groups(P)
	
	for (g, s) in zip(G, S) :
		
		plt.figure(figsize=(9, 3), dpi=300, frameon=False)
		
		# Samplewise classification
		
		m = P[s].as_matrix()
		m.sort(axis=1)
		im = plt.imshow(m.T, aspect='auto', cmap=plt.cm.Blues, vmin=0, vmax=1)
		
		plt.xticks(range(len(P.index)), P.index)
		plt.yticks([], [])
		plt.title(g + " ({} single cells)".format(len(s)))
		
		plt.tight_layout()
		
		# Colorbar
		
		# https://matplotlib.org/users/tight_layout_guide.html
		from mpl_toolkits.axes_grid1 import make_axes_locatable
		cax = make_axes_locatable(plt.gca()).append_axes("right", "2%", pad="3%")
		cb = plt.colorbar(im, cax=cax, ticks=[0, 0.5, 1])
		#
		cb.set_clim(0, 1)
		cax.set_yticklabels(["0%", "confidence", "100%"], va='center', rotation=90, size=5)
		cax.yaxis.set_ticks_position('left')
		cax.tick_params(axis='both', which='both', length=0)
		
		# Save to disk
		
		plt.savefig(OFILE['pred-sc-bygroup'].format(group=g, ext="pdf"))
		plt.close()


def plot_singlecell() :
	(B, M, P) = predict_bcxx(singlecell=True)
	plot_singlecell_samples(P)
	plot_singlecell_subtypes(P)


## == BULK PREDICTIONS == ##

def plot_bulk_subtypes(P) :
	
	for g in P.columns :
		
		plt.figure(figsize=(9, 3), dpi=300, frameon=False)
		
		# Tumorwise classification
		
		plt.bar(range(len(P.index)), list(P[g]), width=1, tick_label=P.index)
		
		plt.xlim((-0.5, len(P.index) - 0.5))
		plt.ylim((0, 1))
		
		plt.yticks([], [])
		plt.title(g + " (bulk)")
		
		plt.tight_layout()
		
		# Colorbar
		
		# https://matplotlib.org/users/tight_layout_guide.html
		from mpl_toolkits.axes_grid1 import make_axes_locatable
		cax = make_axes_locatable(plt.gca()).append_axes("right", "0.01%", pad="5%")
		plt.xticks([], [])
		plt.yticks([0, 0.5, 1], ["0%", "50%", "100%"], rotation=90, size=8)
		
		# Save to disk
		
		plt.savefig(OFILE['pred-bulk'].format(group=g, ext="pdf"))
		plt.close()

def plot_bulk() :
	(B, M, P) = predict_bcxx(singlecell=False)
	plot_bulk_subtypes(P)


## ==================== ENTRY :

if (__name__ == "__main__") :
	plot_sc2bulk()
	plot_singlecell()
	plot_bulk()
	
	# https://github.com/tensorflow/tensorflow/issues/3388#issuecomment-271107725
	keras_backend.clear_session()
