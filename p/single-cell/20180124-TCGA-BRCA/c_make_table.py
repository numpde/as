
# RA, 2018-01-24

# Run as
#    python3 c*.py


## ================== IMPORTS :

import gzip
import sys
import os
import json
import inspect
import pickle
import xml.etree.ElementTree
import pandas as pd

from shutil import copyfile
from glob import glob
from itertools import chain
from collections import defaultdict
from progressbar import ProgressBar as Progress


## ==================== INPUT :

IFILE = {
	'casefile' : "OUTPUT/b_cases/UV/{caseid}/{name}",
}


## =================== OUTPUT :

OFILE = {
	'combined' : "OUTPUT/c_make_table/UV/tcga-brca-fpkm.{ext}",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	# Write the combined table as txt?
	'dump txt' : True,
}

TESTMODE = ("TEST" in sys.argv)


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))


## ===================== WORK :

def get_patient_barcode(clinical) :
	# Get patient barcode from the "clinical" xml file
	r = xml.etree.ElementTree.parse(clinical).getroot()
	patient_barcode = [e for e in r.findall(".//") if e.tag.endswith("patient_barcode")]
	assert(len(patient_barcode) == 1)
	return patient_barcode[0].text


def files_by_patient() :
	
	print("Collecting file paths...")
	
	paths = sorted(glob(IFILE['casefile'].format(caseid="*", name="")))
	
	if TESTMODE : paths = paths[0:3]
	
	P2FF = dict()
	
	for p in Progress()(paths) :
		# Case UUID
		c = os.path.basename(os.path.normpath(p))
		
		def one(L) :
			if not L : return None
		
			assert(type(L) == list)
			assert(len(L) == 1)
			return L[0]
		
		# Files
		F = dict()
		
		F['FPKM']         = glob(p + "*FPKM.txt*")
		F['FPKM-UQ']      = glob(p + "*FPKM-UQ.txt*")
		F['htseq.counts'] = glob(p + "*htseq.counts*")
		
		F['clinical']     = one(glob(p + "*clinical*"))
		F['biospecimen']  = one(glob(p + "*biospecimen*"))
		
		if (not F['clinical']) : 
			print("")
			print("No 'clinical' file in", p)
			continue
		
		patient_barcode = get_patient_barcode(F['clinical'])
		assert(not patient_barcode in P2FF), "Duplicate patient barcode."
		
		P2FF[patient_barcode] = F
		
	return P2FF


def XX_by_patient() :
	
	P2FF = files_by_patient()
	
	print("There are {} patients.".format(len(P2FF)))
	print("No FPKM data found for {} patient(s).".format(sum(not F['FPKM'] for F in P2FF.values())))
	
	P2XX = dict()
	
	print("Collecting FPKM data...")
	
	for (p, F) in Progress()(sorted(P2FF.items())) :
		if not F['FPKM'] : continue
		
		#print("-" * 32)
		#print("Patient:", p)
		#print("Path:", os.path.dirname(F['clinical']))
		
		P2XX[p] = [
			pd.read_table(
				gzip.open(filename, 'rb'), 
				names=['ENSG', "{}_{}".format(p, n)]
			).sort_values('ENSG')
			
			for (n, filename) in enumerate(F['FPKM'])
		]
	
	return P2XX	


def get_combined_table() :
	
	P2XX = XX_by_patient()
	
	print("Combining...")
	
	# Get the set of all genes
	E = set()
	for (p, XX) in P2XX.items() :
		for data in XX :
			E = E | set(data['ENSG'])
	
	combined = pd.DataFrame(data={'ENSG' : sorted(E)})
	
	for (_, XX) in Progress()(sorted(P2XX.items())) :
		for (n, data) in enumerate(XX) :
			combined = pd.merge(combined, data, on='ENSG', how='outer')
	
	return combined


def main() :
	
	X = get_combined_table()
	
	print("Writing files...")
	
	# Write into a pickle file
	pickle.dump(
		{
			'X' : X,
			'script' : THIS,
		},
		open(OFILE['combined'].format(ext="pkl"), 'wb')
	)
	
	# Write into a txt file
	if PARAM['dump txt'] :
		X.to_csv(OFILE['combined'].format(ext="txt"), sep='\t', index=False)


## ==================== ENTRY :

if (__name__ == "__main__") :
	main()
