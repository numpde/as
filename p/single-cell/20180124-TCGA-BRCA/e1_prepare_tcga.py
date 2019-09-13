
# RA, 2018-01-26

# Run as
#    python3 e1*.py PREPARE
#    python3 e1*.py OVERVIEW


## ================== IMPORTS :

import os
import sys
import pickle
import inspect
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## ==================== INPUT :

IFILE = {
	# Expression- and some meta-data
	'TCGA' : "OUTPUT/c_make_table/UV/tcga-brca-fpkm.pkl",
	
	# ENSG to Symbol
	'ENSG' : "OUTPUT/d_ensg/ensg-info.txt",
	
	# Subtype classification
	'subtype' : "ORIGINALS/TCGA-BRCA-01/BRCA.547.PAM50.SigClust.Subtypes.txt",
	
	# TCGA clinical data
	'clinical' : "OUTPUT/b1_clinical/tcga-clinical.pkl",
}


## =================== OUTPUT :

OFILE = {
	# Prepared data
	'DATA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	
	# Compiled classification table
	'class' : "OUTPUT/e_prepared/tcga-class.txt",
	
	# PAM50 classification overview
	'pam50-overview' : "OUTPUT/e_prepared/tcga-pam50-overview/{status}.{ext}",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

PARAM = {
	# Classification methods
	'Classification' : {
		'PAM50' : {
			'key' : 'PAM50',
			'sub' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
		},
		
		'SigClust' : {
			'key' : 'SigClust',
			'sub' : list(range(-13, 1)),
		},
		
		'ER status' : {
			'key' : 'er_status_by_ihc',
			'sub' : ["Positive", "Negative", "Indeterminate"],
		},
		
		'PR status' : {
			'key' : 'pr_status_by_ihc',
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
	},
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Check if pandas series has unique items
def is_unique(S) : return (len(set(S)) == len(S))

# Log which files are opened
def logged_open(filename, mode='r', *argv, **kwargs) :
	print("({}):\t{}".format(mode, filename))
	return open(filename, mode, *argv, **kwargs)


## ========== DATA COLLECTION :

# TCGA expression data
def get_X() :
	
	# Read ENSG info file
	ENSG = pd.read_table(IFILE['ENSG'])
	# Rename columns
	ENSG = ENSG.rename(columns={ 'ensembl_gene_id' : 'ENSG', 'hgnc_symbol' : 'symbol' })
	# Select only those columns
	ENSG = ENSG[['ENSG', 'symbol']]
	# Omit ENSG's without a known gene symbol
	ENSG = ENSG.dropna()
	
	# Read TCGA table
	X: pd.DataFrame
	X = pickle.load(logged_open(IFILE['TCGA'], 'rb'))['X']

	# Drop the version number from the ENSG IDs
	# Do not reindex by ENSG, as other datasets may not have it
	X['ENSG'] = [e[0:15] for e in X['ENSG']]

	# ENSG ID's are unique; gene symbols may not be
	assert(X['ENSG'].is_unique)
	assert(ENSG['ENSG'].index.is_unique)

	# Append the Symbol column, dropping rows where unknown
	X = X.merge(ENSG, on='ENSG', how='inner')

	# At this point, Symbols may repeat, viz (2019-08-03):
	# 	for (symbol, group) in X.groupby('symbol'):
	# 		if (len(group) > 1):
	# 			print(symbol, ":", list(group['ENSG']))
	# yields
	# 	DUXAP8 : ['ENSG00000206195', 'ENSG00000271672']
	# 	GOLGA8M : ['ENSG00000188626', 'ENSG00000261480']
	# 	ITFG2-AS1 : ['ENSG00000256150', 'ENSG00000258325']
	# 	LINC01238 : ['ENSG00000237940', 'ENSG00000261186']
	# 	LINC01505 : ['ENSG00000234229', 'ENSG00000234323']
	# 	PINX1 : ['ENSG00000254093', 'ENSG00000258724']
	# 	POLR2J4 : ['ENSG00000214783', 'ENSG00000272655']
	# 	PRAMEF7 : ['ENSG00000204510', 'ENSG00000279195']
	# 	RMRP : ['ENSG00000269900', 'ENSG00000277027']
	# 	RN7SKP260 : ['ENSG00000252716', 'ENSG00000281402']
	# 	RNU6-771P : ['ENSG00000202250', 'ENSG00000207511']
	# 	SCARNA4 : ['ENSG00000280466', 'ENSG00000281394']
	# 	SNORA16A : ['ENSG00000274582', 'ENSG00000280498']
	# 	SNORA17B : ['ENSG00000276161', 'ENSG00000280496']
	# 	SNORA50A : ['ENSG00000206952', 'ENSG00000281910']
	# 	SNORD38B : ['ENSG00000207421', 'ENSG00000281859']
	# 	SNORD3D : ['ENSG00000277947', 'ENSG00000281000']

	# Index by Symbol, summing over subgroups
	# The ENSG ID column is dropped as a nuisance column:
	# https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
	X = X.groupby('symbol').sum()

	return X


# TCGA subtype classification
def get_class_S() :
	
	# Read the subtype classification
	S = pd.read_table(IFILE['subtype'])
	
	# Rename columns
	S = S.rename(columns={ 
		# Example value
		
		# TCGA-AN-A0FL-01A-11R-A034-07
		'Sample' : 'aliquot_barcode',
		
		# Tumor
		'Type' : 'tissue_type', 
		
		# -13 
		'Siglust' : 'SigClust', # SIC!
		
		# LumA / LumB / Her2 / Basal / Normal
		'PAM50' : 'PAM50',
	})
	
	# Append a column with patient ID (also the case ID)
	S['patient_barcode'] = [a[0:12] for a in S['aliquot_barcode']]
	
	return S


# TCGA clinical information
def get_class_C() :
	
	return pickle.load(logged_open(IFILE['clinical'], 'rb'))['patient_barcode']


def get_class() :
	
	# Sample type for each aliquot
	T = pickle.load(logged_open(IFILE['TCGA'], 'rb'))['FI']
	# Keep only the FPKM file related info
	T = T.loc[ T['type'] == "HTSeq - FPKM", : ]
	# Keep only those columns
	T = T[['aliquot_barcode', 'sample_type']]
	# Sample types
	T['sample_type'] = pd.Categorical(T['sample_type'], ["Primary Tumor", "Solid Tissue Normal", "Metastatic"], ordered=True)
	# Identify record by aliquot barcode
	assert(is_unique(T['aliquot_barcode']))
	# Use the aliquot barcode as index
	T = T.set_index('aliquot_barcode')
	
	
	# Subtypes
	S = get_class_S()
	# Identify record by aliquot barcode
	assert(is_unique(S['aliquot_barcode']))
	# Use the aliquot barcode as index
	S = S.set_index('aliquot_barcode')
	# Some patients have multiple aliquots though
	assert(not is_unique(S['patient_barcode']))
	# But, there is exactly one 'tumor' record per patient
	assert(all((len(S[(S['patient_barcode'] == p) & (S['tissue_type'] == 'tumor')]) == 1) for p in S['patient_barcode']))
	
	
	# Clinical
	C = get_class_C()
	# Identify record by patient barcode
	assert(is_unique(C.columns))
	
	# Classes info
	CL = PARAM['Classification']
	
	# Merged classification table
	M = pd.DataFrame(index=T.index, columns=CL.keys())
	
	# Make columns categorical
	for c in M.columns :
		M[c] = pd.Categorical(M[c], categories=CL[c]['sub'], ordered=True)
	
	# Append tissue type column
	M['sample_type'] = T
	
	# Collect available classification for each aliquot
	for a in M.index :
		
		# Fetch PAM50 / SigClust data
		if (a in S.index) :
			for c in M.columns :
				if (c in CL) :
					k = CL[c]['key']
					if (k in S.columns) :
						M.loc[a, c] = S[k][a]
		
		# Skip aliquot if we know it is not from (primary) tumor ...
		if (a in T.index) and ("umor" not in T['sample_type'][a]) : continue
		if (a in S.index) and ("umor" not in S['tissue_type'][a]) : continue
		
		# ... otherwise extract clinical info by *patient* barcode
		
		# Patient barcode
		p = a[0:12]
		
		# Have clinical data?
		if (p in C.columns) : 
			for c in M.columns :
				if (c in CL) : 
					k = CL[c]['key']
					if (k in C.index) :
						M.loc[a, c] = C[p][k]
	
	M.to_csv(OFILE['class'], sep='\t')
	
	return M


def PREPARE() :
	
	pickle.dump(
		{
			# Expression table
			'X' : get_X(),
			
			# For backwards compatibility (2018-02-04)
			'subtype' : get_class_S(), 
			
			# Classification data
			'class' : get_class(),
			'class-meta' : PARAM['Classification'],
			
			# This file
			'script' : THIS,
		},
		logged_open(OFILE['DATA'], 'wb')
	)


## ================= PLOTTING :

def plot_pam50_overview(M) :
	# Compare PAM50 and ER/HER
	
	C = ['ER status', 'PR status', 'HER2 status (IHC)', 'HER2 status (FISH)', 'TN']
	
	def short(t) :
		if (type(t) is float) : return t
		if (type(t) is str) :
			t = t.replace("Positive", "(+)")
			t = t.replace("Negative", "(--)")
			t = t.replace("Equivocal", "(~)")
			t = t.replace("Indeterminate", "(?)")
			return t
		else :
			return [short(s) for s in t]
	
	# Make local copy for modification
	M = M.copy()
	
	# Add a triple-negative status column
	M['TN'] = pd.Categorical(["N/A"] * M.index.size, categories=["Yes", "No", "Normal", "N/A"], ordered=True)
	#
	for (i, r) in M.iterrows() :
		if ("ormal" in r['sample_type']) :
			M.loc[i, 'TN'] = "Normal"
			continue
			
		if (r['ER status'] == r['PR status'] == r['HER2 status (IHC)'] == "Negative") :
			M.loc[i, 'TN'] = "Yes"
			continue
		
		if ("Positive" in [r['ER status'], r['PR status'], r['HER2 status (IHC)']]) :
			M.loc[i, 'TN'] = "No"
			continue
	
	for (n, c) in enumerate(C) :
		# plt.figure()

		fig: plt.Figure
		ax: plt.Axes
		(fig, ax) = plt.subplots(figsize=(4, 3), dpi=100)

		p = 'PAM50'
		m = M[[p, c]].dropna()
		S = pd.DataFrame(0, index=M[p].cat.categories, columns=M[c].cat.categories)
		for (_, r) in m.iterrows() : S.loc[r[p], r[c]] += 1
		
		def norm(s) : return s / (sum(s) or 1)
		
		im = ax.imshow(S.apply(norm, axis=0), cmap=plt.cm.Blues, aspect='auto', vmin=0, vmax=1)
		
		for (x, j) in enumerate(S.columns) :
			for (y, i) in enumerate(S.index) :
				ax.text(x, y, S.loc[i, j], va='center', ha='center')

		ax.tick_params(axis='both', which='both', length=0)

		(xlim, ylim) = (ax.get_xlim(), ax.get_ylim())

		ax.set_xticks(range(0, len(S.columns)))
		ax.set_xticklabels(short(S.columns))

		ax.set_yticks(range(0, len(S.index)))
		ax.set_yticklabels(S.index)

		ax.set_xlabel(c)
		ax.set_ylabel("Classified (in TCGA) as...")

		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
	
		fig.tight_layout()
		
		# Colorbar

		# https://stackoverflow.com/questions/13784201/matplotlib-2-subplots-1-colorbar
		(y0, y1) = (ax.get_position().y0, ax.get_position().y1)
		cax = fig.add_axes([0.98, y0, 0.01, y1-y0])
		fig.subplots_adjust(right=0.9)
		cb = fig.colorbar(im, cax=cax, ticks=[0, 0.5, 1], )
		cax.tick_params(labelsize=5) 
		im.set_clim(0, 1)
		cb.set_ticks([0, 0.5, 1])
		cax.set_yticklabels(["0%", "proportion", "100%"], va='center', rotation=90)
		cax.yaxis.set_ticks_position('left')
		cax.tick_params(axis='both', which='both', length=0)
		
		# Save to disk

		for ext in ["png", "pdf"]:
			fig.savefig(
				OFILE['pam50-overview'].format(status=(c.replace(" ", "-")), ext=ext),
				bbox_inches='tight', pad_inches=0,
				dpi=300
			)
		plt.close(fig)


def OVERVIEW() :
	DATA = pickle.load(open(OFILE['DATA'], 'rb'))
	plot_pam50_overview(DATA['class'])

## ==================== ENTRY :

if (__name__ == "__main__") :
	if ("PREPARE"  in sys.argv) : PREPARE()
	if ("OVERVIEW" in sys.argv) : OVERVIEW()

