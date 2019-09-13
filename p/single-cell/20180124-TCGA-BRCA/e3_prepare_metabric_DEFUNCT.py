
# RA, 2018-02-06

# Run as
#    python3 e3*.py


## ================== IMPORTS :

import os
import re
import pickle
import inspect
import tarfile, io, urllib.request
import pandas as pd


## ==================== INPUT :

IFILE = {
	# METABRIC data, if available
	'BRIC' : "ORIGINALS/METABRIC/UV/brca_metabric.tar.gz",
	
	# Otherwise, download from here
	'BRIC URL' : "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric.tar.gz",
}


## =================== OUTPUT :

OFILE = {
	# Put the downloaded file here, if not available
	'BRIC' : IFILE['BRIC'],
	
	# Prepared data
	'DATA' : "OUTPUT/e_prepared/UV/bric.pkl",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

PARAM = {
	
	# Files in the archive
	'data-expression' : "brca_metabric/data_expression.txt",
	'meta-expression' : "brca_metabric/meta_expression.txt",
	'data-clinical'   : "brca_metabric/data_clinical_patient.txt",
	'meta-clinical'   : "brca_metabric/meta_clinical_patient.txt",
	
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


## ===================== WORK :

def prepare() :
	
	# Obtain the tarball if necessary
	if not os.path.isfile(IFILE['BRIC']) :
		
		with urllib.request.urlopen(IFILE['BRIC URL']) as response :
			with logged_open(OFILE['BRIC'], 'wb') as f :
				f.write(response.read())
		
		assert(os.path.isfile(IFILE['BRIC']))
	
	# Extract data from tarball
	with tarfile.open(fileobj=logged_open(IFILE['BRIC'], 'rb'), mode='r:gz') as t : 
		
		with t.extractfile(t.getmember(PARAM['data-expression'])) as m :
			X = pd.read_csv(io.BytesIO(m.read()), sep='\t')
		
		with t.extractfile(t.getmember(PARAM['meta-expression'])) as m :
			x = m.read()
			
		with t.extractfile(t.getmember(PARAM['data-clinical'])) as m :
			C = pd.read_csv(io.BytesIO(m.read()), sep='\t', comment='#')
		
		with t.extractfile(t.getmember(PARAM['meta-clinical'])) as m :
			c = m.read()
	
	# Rename column
	X = X.rename(columns={ "Hugo_Symbol" : "Symbol" })
	# Keep only the single sample columns
	X = X[['Symbol'] + sorted(c for c in X if re.match("MB-[0-9]+", c))]
	# Reindex by Symbol, summing over subgroups
	assert(is_unique(X['Symbol']))
	X = X.set_index('Symbol')
	
	# Save to disk
	
	pickle.dump(
		{
			'X' : X,
			'C' : C,
			'X-meta' : x,
			'C-meta' : c,
			'script' : THIS,
		},
		logged_open(OFILE['DATA'], 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	prepare()


