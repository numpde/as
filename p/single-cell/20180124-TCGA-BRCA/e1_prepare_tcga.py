
# RA, 2018-01-26

# Run as
#    python3 e1*.py


## ================== IMPORTS :

import os
import pickle
import inspect
import pandas as pd


## ==================== INPUT :

IFILE = {
	# Extract the list of relevant genes from here
	'TCGA' : "OUTPUT/c_make_table/UV/tcga-brca-fpkm.pkl",
	
	# ENSG to Symbol
	'ENSG' : "OUTPUT/d_ensg/ensg-info.txt",
	
	# Subtype classification
	'subtype' : "ORIGINALS/TCGA-BRCA-01/BRCA.547.PAM50.SigClust.Subtypes.txt",
}


## =================== OUTPUT :

OFILE = {
	# Prepared data
	'DATA' : "OUTPUT/e_prepared/UV/tcga.pkl",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Log which files are written
def logged_open(filename, mode='r', *argv, **kwargs) :
	print("({}):\t{}".format(mode, filename))
	return open(filename, mode, *argv, **kwargs)

## ===================== WORK :

def prepare() :
	
	# Read ENSG info file
	ENSG = pd.read_table(IFILE['ENSG'])
	# Rename columns
	ENSG = ENSG.rename(columns={ 'ensembl_gene_id' : 'ENSG', 'hgnc_symbol' : 'symbol' })
	# Select only those columns
	ENSG = ENSG[['ENSG', 'symbol']]
	# Omit rows with missing data
	ENSG = ENSG.dropna(axis=0, how='any')

	# Read TCGA table
	X = pickle.load(open(IFILE['TCGA'], 'rb'))['X']
	# Omit version number from the ENSG IDs
	X['ENSG'] = [e[0:15] for e in X['ENSG']]
	# Append the Symbol column, dropping rows where unknown
	# Then "ENSG" is the index column
	X = X.merge(ENSG, left_index=True, right_index=True, how='inner')
	# Index by Symbol, summing over subgroups
	# The ENSG index column is dropped
	X = X.groupby('symbol').sum()
	
	# Read the subtype classification
	C = pd.read_table(IFILE['subtype'])
	# Rename columns
	C = C.rename(columns={ 
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
	C['patient_barcode'] = [s[0:12] for s in C['aliquot_barcode']]
	
	pickle.dump(
		{
			'X' : X,
			'subtype' : C,
			'script' : THIS,
		},
		logged_open(OFILE['DATA'], 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	prepare()

