
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
	
	'clinical' : "OUTPUT/b1_clinical/tcga-clinical.pkl",
}


## =================== OUTPUT :

OFILE = {
	'classified' : "OUTPUT/h2_class/TCGA-BCXX/{method}/{geneset}.{ext}",
}


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Clinical classification methods
	'Classification' : {
		'ER status' : {
			'key' : 'er_status_by_ihc',
			'sub' : ["Positive", "Negative", "Indeterminate"],
		},
		
		'HER2 status (FISH)' : {
			'key' : 'her2_fish_status',
			'sub' : ["Positive", "Negative", "Equivocal", "Indeterminate"],
		},
		
		'HER2 status (IHC)' : {
			'key' : 'her2_status_by_ihc',
			'sub' : ["Positive", "Negative", "Equivocal", "Indeterminate"],
		},
	}
}


## ====================== AUX :

# Change GO:000 to GO-000 in filenames
def nicer(filename) : return filename.replace("GO:", "GO-")

def abbr(L) :
	if (type(L) is str) : 
		ab = {
			"Positive" : "(+)",
			"Negative" : "(--)",
			"Equivocal" : "(~)",
			"Indeterminate" : "(?)",
		}
		for (a, b) in ab.items() : L = L.replace(a, b)
		return L
	else :
		# Guess that L is an iterable
		return [abbr(i) for i in L]

## ====================== (!) :

pass


## ===================== WORK :

# Compute the similarity grouping by class
def classify(M, method) :
	
	# Clinical classification key in datatable
	clinikey = PARAM['Classification'][method]['key']
	# Possible subtypes
	subtypes = PARAM['Classification'][method]['sub']
	
	# Note:
	# M is a (TCGA bulk sample) x (BCXX single cell sample) similarity matrix
	
	# Partial clinical classification of TCGA samples
	C = pickle.load(open(IFILE['clinical'], 'rb'))['patient_barcode']
	# Select the field of the clinical classification
	C = C.transpose()[clinikey].dropna()
	# Classification subtypes are represented
	assert(set(subtypes) == set(C.unique()))
	
	# Aliquot to class
	def a2c(a) : 
		a = a[0:len(C.index[0])]
		if a in C : return C[a]
		return np.nan
	
	# 
	C = pd.Series(M.index, index=M.index).apply(a2c).dropna()

	# Keep only the TCGA samples that have been classified
	M = M.loc[ list(pd.Series(M.index).isin(C.index)), : ]
	
	# Group into a (TCGA subtype) x (BCXX tumor) matrix
	def group(M) :
		
		# Omit records with missing data
		N = M.dropna().astype(float)
		
		# Group by class
		N = N.groupby(
			# Aliquot barcode to class
			lambda a : C[ a ],
			axis=0
		).mean()
		
		# Group by BCXX tumors
		N = N.groupby(
			# Sample "BCXX[LN][_Re]_YY" to cluster "BCXX"
			lambda c : c[:-3],
			axis=1
		).mean()
		
		# Use the order of the subtypes list
		N = N.reindex(pd.Categorical(N.index, categories=subtypes)).sort_index()
		
		return N
	
	# Group into a (TCGA subtype) x (BCXX tumor) matrix
	N = group(M)

	return N

# Plot a (TCGA type) x (BCXX tumor) matrix
def plot(N):

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
	ax.set_yticklabels(abbr(N.index))
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


def COMPUTE() :
	
	all_cossim = pickle.load(open(IFILE['cossim'], 'rb'))['cossim']
	
	for cossim in commons.Progress()(all_cossim) :
		for method in PARAM['Classification'].keys() :
			meta = cossim['meta']
			setid = meta['id']
			
			#print("Gene set:", setid)
			#print("Method:", method)

			try:
				N = classify(cossim['M'], method)
			except:
				commons.logger.warning("Failed on gene set:", setid)

			(fig, ax) = plot(N)

			title = "{name} ({len} genes)".format(name=meta['info'], len=len(meta['set']))
			ax.set_title(title, fontsize=6)

			ax.set_ylabel(method + " normed cos-sim to...")

			for ext in ['png', 'pdf']:
				fig.savefig(
					commons.makedirs(OFILE['classified'].format(method=method, geneset=nicer(setid), ext=ext)),
					bbox_inches='tight', pad_inches=0,
					dpi=180
				)

			plt.close(fig)


## ==================== ENTRY :

if (__name__ == "__main__") :
	COMPUTE()

