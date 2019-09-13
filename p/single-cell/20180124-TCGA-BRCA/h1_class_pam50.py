
# RA, 2018-01-28


## ================== IMPORTS :

from helpers import commons

import os
import sys
import math
import pickle
import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## ==================== INPUT :

IFILE = {
	'cossim' : "OUTPUT/g_cossim/UV/TCGA-v-BCXX.pkl",
	
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'clinical' : None,
}


## =================== OUTPUT :

OFILE = {
	'classified' : "OUTPUT/h1_class/TCGA-BCXX/PAM50/{geneset}.{ext}",
}


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# PAM50 label order
	'PAM50' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
}


## ====================== AUX :

# Change GO:000 to GO-000 in filenames
def nicer(filename) : return filename.replace("GO:", "GO-")

## ====================== (!) :

pass


## =================== WORKER :

# Compute the similarity grouping by class
def classify(M) :
	
	# Note:
	# M is a (TCGA bulk sample) x (BCXX single cell sample) similarity matrix
	
	# Partial classification of TCGA samples, from
	# http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt
	C = pickle.load(open(IFILE['TCGA'], 'rb'))['subtype']
	# These are the PAM50 subtypes
	PAM50 = PARAM['PAM50']
	# Assume uniqueness of 'aliquot_barcode's for unambiguous classification
	assert(C['aliquot_barcode'].is_unique)
	# Expect all the PAM50 subtypes to be represented
	assert(set(PAM50) == set(C['PAM50']))

	# Index by aliquot_barcode
	C = C.set_index('aliquot_barcode')
	
	# Keep only the TCGA samples that have been PAM50-classified
	M = M.loc[ M.index.isin(C.index), : ]
	# Each sample is now classified
	assert(set(C.index) >= set(M.index))

	# Group into a (TCGA PAM50 type) x (BCXX tumor) matrix
	def group(M) :
		
		# Omit records with missing data
		N = M.dropna().astype(float)

		# Group by PAM50 type
		N = N.groupby(
			# Aliquot barcode to PAM50 type
			lambda a : C.loc[a, 'PAM50'],
			axis=0
		).mean()

		# Group by BCXX tumors
		N = N.groupby(
			# Sample "BCXX[LN][_Re]_YY" to cluster "BCXX"
			lambda c : c[:-3],
			axis=1
		).mean()

		# Use the order of the PAM50 list
		N = N.reindex(pd.Categorical(N.index, categories=PAM50)).sort_index()
		
		return N

	# Group into a (TCGA PAM50 type) x (BCXX tumor) matrix
	N = group(M)

	return N


# Plot a (TCGA type) x (BCXX tumor) matrix
def plot(N):

	N = N.drop(index='Normal')

	# Normalize column-wise
	for c in N.columns : N[c] /= N[c].sum()

	## Normalize row-wise
	#for i in N.index : N.ix[i] /= N.ix[i].sum()

	fig: plt.Figure
	ax: plt.Axes
	(fig, ax) = plt.subplots()

	im = ax.imshow(N.values, aspect='auto', cmap=plt.cm.Blues)
	(xlim, ylim) = (ax.get_xlim(), ax.get_ylim())

	ax.set_xticks(list(range(len(N.columns))))
	ax.set_xticklabels(N.columns, rotation=60)
	ax.set_yticks(list(range(len(N.index))))
	ax.set_yticklabels(N.index)
	ax.tick_params(axis='x', which='major', labelsize=6)
	ax.tick_params(axis='y', which='major', labelsize=8)
	ax.set_ylabel("Normed cos-sim to...")

	ax.set_xlim(xlim)
	ax.set_ylim(ylim)

	fig.tight_layout()

	# https://stackoverflow.com/questions/13784201/matplotlib-2-subplots-1-colorbar
	(y0, y1) = (ax.get_position().y0, ax.get_position().y1)
	cax = fig.add_axes([0.98, y0, 0.01, y1 - y0])
	fig.subplots_adjust(right=0.9)
	cb = fig.colorbar(im, cax=cax, ticks=[0, 0.5, 1], )
	cax.tick_params(labelsize=5)
	im.set_clim(0, 1)
	cb.set_ticks([0, 0.5, 1])
	cax.set_yticklabels(["0%", "proportion", "100%"], va='center', rotation=90)
	cax.yaxis.set_ticks_position('left')
	cax.tick_params(axis='both', which='both', length=0)

	return (fig, ax)


## =================== MASTER :

def main() :
	
	all_cossim = pickle.load(open(IFILE['cossim'], 'rb'))['cossim']
	
	for cossim in commons.Progress()(all_cossim) :
		setid = cossim['meta']['id']

		try:
			N = classify(cossim['M'])
		except:
			commons.logger.warning("Failed on gene set: {}".format(setid))
			continue

		(fig, ax) = plot(N)
		title = "{name} ({len} genes)".format(name=cossim['meta']['info'], len=len(cossim['meta']['set']))
		ax.set_title(title, fontsize=6)

		for ext in ['png', 'pdf']:
			fig.savefig(
				commons.makedirs(OFILE['classified'].format(geneset=nicer(setid), ext=ext)),
				bbox_inches='tight', pad_inches=0,
				dpi=180
			)

		plt.close(fig)


## ==================== ENTRY :

if (__name__ == "__main__") :
	main()

