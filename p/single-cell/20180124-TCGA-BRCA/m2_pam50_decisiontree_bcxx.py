
# RA, 2019-08-02

# Apply the mock PAM50 classifier to the single-cell BCXX dataset


from helpers import commons

import pandas as pd
import numpy as np
import json, pickle

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from itertools import chain

from sklearn.ensemble import RandomForestClassifier


PARAM = {
	# Input: Preprocessed BCXX gene expression dataset (GSE75688)
	'BCXX': "OUTPUT/e_prepared/UV/bcxx.pkl",

	# Input: Regression model
	'model': "OUTPUT/m1_pam50_decisiontree/rf_model.dat",

	# List of PAM50 genes
	'PAM50': json.load(open("ORIGINALS/GeneSets/pam50/pam50.json", 'r'))['data'],

	# PAM50 labels order and plot markers
	'PAM50-labels' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
	'PAM50-markers': ['o', '^', 'v', '*', 's'],

	# Normalize the fancy plot to "max = 1"
	'normalize-y': True,

	# Output: Barplot tumor-wise as predicted by classifier
	'pred_classified': "OUTPUT/m2_pam50_decisiontree_bcxx/hist_{kind}.{ext}",
	# Output: Confidences for each cell / each class
	'full_classified': "OUTPUT/m2_pam50_decisiontree_bcxx/full_hist.{ext}",
}


def plot_fancy(Y: pd.DataFrame):

	# Maps "BC06_12" to "BC06"
	parent = (lambda sample_name: sample_name[0:-3])
	# List of tumor labels
	tumors = sorted(set(map(parent, Y.index)))

	# PAM50 types
	classes = Y.columns

	# Class --> color
	colors = dict(zip(PARAM['PAM50-labels'], map(mcolors.to_rgb, mcolors.TABLEAU_COLORS.values())))

	# Sample label proba to (predicted class, -confidence)
	def class_confidence(p):
		n = np.argmax(p.values)
		return (n, -p[n])

	# Max-normalize the whole matrix
	if PARAM['normalize-y']:
		commons.logger.info("Max Y: {}".format(max(Y.values.flatten())))
		Y /= max(Y.values.flatten())

	# Split Y by tumor
	TY = [Y[Y.index.map(parent) == t] for t in tumors]
	del Y

	# Sort sample groups by predicted class and confidence
	TY = [Y.reindex(sorted(Y.index, key=(lambda i: class_confidence(Y.loc[i])))) for Y in TY]

	# Size of each cluster
	R = [len(Y.index) for Y in TY]
	# Insert spacing between them
	spacing = 2
	R = list(chain(*zip([spacing] * len(R), R)))[1:]
	# Add space for the subtypes labels (and the colorbar)
	# R = [sum(R) * 0.03] + [spacing] + R + [sum(R) * 0.03, sum(R) * 0.01, spacing] # with colorbar
	R = [sum(R) * 0.03] + [spacing] + R + [spacing, sum(R) * 0.01, spacing]
	R = np.cumsum(R)
	R = R / max(R)

	with plt.style.context("dark_background"):
		fig = plt.figure(figsize=(10, 2))

		# Vertical orientation
		origin = ["upper", "lower"][1]

		# Ploting axes for the groups
		AX = [
			fig.add_axes([a, 0, b - a, 1], alpha=1)
			for (a, b, _) in zip(R[1::2], R[2::2], tumors)
			# The third component is to have the correct number of axes
		]

		# Axes frame edge color
		for ax in AX:
			plt.setp(ax.spines.values(), color='C6')

		# Group-wise sample classification

		def to_rgba(Y: pd.DataFrame):
			C = [
				list(map((lambda a: [*colors[pam], a]), Y[pam].values))
				for pam in Y.columns
			]
			return np.asarray(C)

		for (t, ax, Y) in zip(tumors, AX, TY):

			im = ax.imshow(
				to_rgba(Y),
				aspect='auto', origin=origin,
				cmap=plt.cm.get_cmap('Blues'),
				vmin=0, vmax=1
			)

			ax.set_xticks([])
			ax.set_yticks([])

			# Group label
			ax.text(
				0.95, 0.05, t,
				size=6,
				transform=ax.transAxes, rotation=90,
				ha='right', va='bottom',
				color="white"
			)

		# # Colorbar
		#
		# cax = fig.add_axes([R[-3], 0.05, R[-2] - R[-3], 0.9], frame_on=False)
		# cb = plt.colorbar(im, cax=cax, ticks=[0, 0.5, 1])
		# #
		# im.set_clim(0, 1)
		# cax.set_yticklabels(["0%", "confidence", "100%"], va='center', rotation=90, size=5)
		# cax.yaxis.set_ticks_position('left')
		# cax.tick_params(axis='both', which='both', length=0)

		# Class labels

		ax = fig.add_axes([0, 0, R[0], 1])
		ax.axis('off')
		#
		L = list(classes)
		if (origin == "upper"): L = list(reversed(L))
		#
		for (y, t) in zip(np.linspace(0, 1, 1 + 2 * len(L))[1::2], L):
			ax.text(
				0.5, y, t,
				size=8,
				transform=ax.transAxes, rotation=90,
				ha='center', va='center',
			)

		# Save figure

		for ext in ['png', 'pdf']:
			fig.savefig(
				commons.makedirs(PARAM['full_classified'].format(ext=ext)),
				transparent=False,
				facecolor=(0.1, 0.1, 0.1, 1),
				# frameon=True,
				bbox_inches='tight', pad_inches=0.1,
				dpi=300
			)

		plt.close(fig)


def plot_hist(y: pd.Series):

	# Maps "BC06_12" to "BC06"
	parent = (lambda sample_name: sample_name[0:-3])
	# List of tumor labels
	tumors = sorted(set(map(parent, y.index)))

	# Histogram (Tumor) x (PAM50 type)
	h = pd.DataFrame({'Tumor': map(parent, y.index), 'Type': y.values}).pivot_table(index='Tumor', columns='Type', aggfunc='size', fill_value=0)

	colors = {pam: ("C{}".format(n)) for (n, pam) in enumerate(PARAM['PAM50-labels'])}

	for kind in ['abs', 'rel']:

		if (kind == 'rel'):
			h = h.apply(lambda c: 100 * c / c.sum(), axis=1)
			zorder = 0
		else:
			zorder = 10

		fig: plt.Figure
		ax: plt.Axes
		(fig, ax) = plt.subplots()

		bottom = 0
		for pam in reversed(PARAM['PAM50-labels']):
			ax.bar(tumors, h[pam], bottom=bottom, label=pam, zorder=zorder, color=colors[pam])
			bottom += h[pam]

		if (kind == 'abs'):
			# zorder doesn't seem to help here
			ax.grid(axis='y', zorder=-10)

		ax.tick_params(axis='x', which='both', labelsize='small', length=0, rotation=45)

		if (kind == 'rel'):
			ax.tick_params(labelleft=False)
			# ax.set_yticklabels(ax.get_yticklabels())


		# Legend in inverse order
		ax.legend(*map(reversed, ax.get_legend_handles_labels()), loc='upper right')

		fig.savefig(
			commons.makedirs(PARAM['pred_classified'].format(kind=kind, ext='png')),
			bbox_inches='tight', pad_inches=0,
			dpi=300
		)

		plt.close(fig)


def classify_and_plot():
	# Expression table (HGNC Symbol x Sample ID)
	# Transpose to a table (Sample ID x HGNC Symbol)
	# Normalize to z-score sample-wise
	# Restrict features to only the PAM50 genes
	X: pd.DataFrame
	X = commons.zscore(pickle.load(open(PARAM['BCXX'], 'rb'))['X'].T, axis=1)[PARAM['PAM50']]

	# Classifier
	model: RandomForestClassifier
	model = pickle.load(open(PARAM['model'], 'rb'))['model']

	# Check that the data layout is compatible with the classifier
	assert(all(list(X.columns == pickle.load(open(PARAM['model'], 'rb'))['features'])))

	# Prediction confidence table (Sample) x (Label)
	Y = pd.DataFrame(data=model.predict_proba(X), index=X.index, columns=model.classes_)
	plot_fancy(Y)

	# Predicted label for each sample in X
	y = pd.Series(model.predict(X), index=X.index)
	plot_hist(y)


if (__name__ == "__main__"):
	classify_and_plot()
