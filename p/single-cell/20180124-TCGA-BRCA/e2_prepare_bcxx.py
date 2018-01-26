
# RA, 2018-01-26

# Run as
#    python3 e2*.py


## ================== IMPORTS :

import os
import re
import pickle
import inspect
import gzip, urllib.request
import pandas as pd


## ==================== INPUT :

IFILE = {
	# BCXX data, if available
	'BCXX' : "ORIGINALS/GSE75688-BCXX/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",
	
	# Otherwise, download from here
	'BCXX URL' : "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75688/suppl/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz",
}


## =================== OUTPUT :

OFILE = {
	# Put the downloaded file here, if not available
	'BCXX' : IFILE['BCXX'],
	
	# Prepared data
	'DATA' : "OUTPUT/e_prepared/UV/bcxx.pkl",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))


## ===================== WORK :

def prepare() :
	
	if not os.path.isfile(IFILE['BCXX']) :
		with urllib.request.urlopen(IFILE['BCXX URL']) as response :
			with open(OFILE['BCXX'], 'wb') as f :
				f.write(gzip.decompress(response.read()))
	
	# Read BC DATA
	X = pd.read_table(IFILE['BCXX'])
	# Rename column
	X = X.rename(columns={ "gene_name" : "Symbol" })
	# Keep only those columns
	X = X[["Symbol"] + [c for c in X if re.match("BC[0-9]+_[0-9]+", c)]]
	# Reindex by Symbol, summing over subgroups
	X = X.groupby("Symbol").sum()
	
	
	pickle.dump(
		{
			'X' : X,
			'script' : THIS,
		},
		open(OFILE['DATA'], 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	prepare()


