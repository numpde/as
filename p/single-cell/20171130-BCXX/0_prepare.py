
# RA, 2019-02-12

# Run as
#    python3 0_prepare.py


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
	'BCXX' : "ORIGINALS/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",
	
	# Otherwise, download from here
	'BCXX URL' : "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75688/suppl/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz",
	
	# Convert possible synonyms to the primary HGNC symbol
	'Syn=>Symb' : "ORIGINALS/go-1/syn2symb.txt",
}


## =================== OUTPUT :

OFILE = {
	# Put the downloaded file here, if not available
	'BCXX' : IFILE['BCXX'],
	
	# Prepared data
	'DATA' : "OUTPUT/0_prepare/UV/bcxx.pkl"
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
	# Keep only the single sample columns
	X = X[['gene_name'] + [c for c in X if re.match("BC[0-9]+.*_[0-9]+", c)]]
	
	# Append a "Symbol" column, converting possible synonyms to pripary HGNC symbols
	#
	# Syn : Symbol synonym --> Primary symbol
	Syn = {
		r[0] : r[1]
		for r in pd.read_csv(IFILE['Syn=>Symb'], sep='\t', index_col=0).itertuples()
	}
	#
	X['Symbol'] = X['gene_name']
	#
	for i in (set(X.index) & set(Syn.keys())) : 
		X.loc[i, 'Symbol'] = Syn[i]
	
	# Reindex by Symbol, summing over subgroups
	X = X.groupby('Symbol').sum()
	
	# Group samples by patient BCXX
	P2SH = {
		p : [(s, h) for (s, h) in enumerate(X.columns) if re.match(p + "_[0-9]+", h)]
		for p in [h[0:4] for h in X.columns]
	}

	# Groups samples by batch BCXX[LN][_Re]
	B2SH = {
		b : [(s, h) for (s, h) in enumerate(X.columns) if re.match(b + "_[0-9]+", h)]
		for b in [h[:-3] for h in X.columns]
	}
	
	# Save to disk
	pickle.dump(
		{
			'X'         : X,    # Gene expression data
			'P2SH'      : P2SH, # Group by patient
			'B2SH'      : B2SH, # Group by batch
			'script'    : THIS, # This script
		},
		open(OFILE['DATA'], 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	prepare()

