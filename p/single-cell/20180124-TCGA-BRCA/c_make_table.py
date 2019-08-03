
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
	'casefile' : "OUTPUT/a_cases/UV/{caseid}/{name}",
	'filemeta' : "ORIGINALS/TCGA-BRCA-01/transcriptome_list.tsv",
}


## =================== OUTPUT :

OFILE = {
	'combined' : "OUTPUT/c_make_table/UV/tcga-brca-fpkm.{ext}",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)

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
		
		# Expects an empty list or a singleton
		def onlyone(L) :
			if not L : return None
		
			assert(type(L) == list)
			assert(len(L) == 1)

			return L[0]
		
		# Files
		F = dict()
		
		F['FPKM']         = glob(p + "*FPKM.txt*")
		F['FPKM-UQ']      = glob(p + "*FPKM-UQ.txt*")
		F['htseq.counts'] = glob(p + "*htseq.counts*")
		
		F['clinical']     = onlyone(glob(p + "*clinical*.xml"))
		F['biospecimen']  = onlyone(glob(p + "*biospecimen*"))
		
		if not F['clinical'] : 
			print("Note: No 'clinical' file in", p)
		
		# Not all cases have the 'clinical' file, so cannot use
		#     get_patient_barcode(F['clinical'])
		# Use the case UUID instead:
		patient_uuid = c
		
		assert(not patient_uuid in P2FF), "Duplicate patient barcode."
		
		P2FF[patient_uuid] = F
		
	return P2FF


def get_transcriptome_file_info() :

	FI = pd.read_table(IFILE['filemeta'])

	FI.rename(inplace=True, columns={
		# [Found in biospecimen / clinical file?] Example(s):
		
		# [Y/N] TCGA-E9-A3QA-01A-61R-A22K-07
		'cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id' : 
			'aliquot_barcode', 
		
		# [Y/N] 8a7e272c-39f2-42e6-b1d5-7225cf12fceb
		'cases.0.samples.0.portions.0.analytes.0.aliquots.0.aliquot_id' : 
			'aliquot_uuid',
		
		# [Y/N] Primary Tumor / Solid Tissue Normal / Metastatic
		'cases.0.samples.0.sample_type' : 
			'sample_type',  
		
		# [Y/N] 27e8b056-063f-4742-a591-b98fb9a948a5
		'cases.0.samples.0.sample_id' : 
			'sample_uuid',  
		
		# [Y/Y] TCGA-E9-A3QA
		'cases.0.submitter_id' : 
			'patient_barcode', 
		
		# [Y/Y] f3cb557d-23e4-4fd1-81ca-db1a3f56d56e
		'cases.0.case_id' : 
			'patient_uuid', 
		
		# [N/N] 148d950b-4202-4f3f-be15-84735cd08a48.FPKM.txt.gz
		'file_name' : 
			'file_name', 
		
		# [N/N] HTSeq - FPKM
		'analysis.workflow_type' : 
			'type', 
		
		# [N/N] Transcriptome Profiling
		'data_category' : 
			'data_category',
		
		# [N/N] Gene ... / Isoform ... / miRNA Expression Quantification
		'data_type' : 
			'data_type', 
		
		# [N/N] 61b737e2-5a35-4c17-b674-f4ddce8355bb
		'file_id' : 
			'file_uuid', 
	})
	
	# Remove the rendundant column 'id'
	assert((FI['id'] == FI['file_uuid']).all())
	del FI['id']

	return FI


def XX_by_patient() :
	
	FI = get_transcriptome_file_info()

	#
	assert(FI['file_name'].is_unique)
	FI = FI.set_index('file_name')
	
	P2FF = files_by_patient()
	
	print("There are {} patients.".format(len(P2FF)))
	print("No FPKM data found for {} patient(s).".format(sum(not F['FPKM'] for F in P2FF.values())))
	
	P2XX = defaultdict(list)
	
	print("Collecting FPKM data...")
	
	for (p, F) in Progress()(sorted(P2FF.items())) :
		if not F['FPKM'] : continue
		
		#print("-" * 32)
		#print("Patient:", p)
		#print("Path:", os.path.dirname(F['clinical']))
		
		for (n, filename) in enumerate(F['FPKM']) :

			i = os.path.basename(filename)

			if i not in FI.index:
				print("(!) File {i} not in FI".format(i=i))
				continue
			
			# Constency of patient UUID info coming from the directory tree and the 'filemeta' list
			patient_uuid = FI.loc[ i, 'patient_uuid' ]
			assert(p == patient_uuid), "Patient ID mismatch"
			
			# Aliquot ID of this file (file of different data kinds that belong together share this ID)
			aliquot_barcode = FI.loc[ i, 'aliquot_barcode' ]
			
			P2XX[p].append(
				pd.read_table(
					# The file contains two columns
					gzip.open(filename, 'rb'), 
					# These will be the headers of the columns
					names = [ 'ENSG', aliquot_barcode ]
				).sort_values('ENSG')
			)
	
	P2XX = dict(P2XX)
	
	return (FI, P2XX)


def get_combined_table() :
	
	(FI, P2XX) = XX_by_patient()
	
	print("Combining...")
	
	# Get the set of all genes
	E = set()
	for (p, XX) in P2XX.items() :
		for data in XX :
			E = E | set(data['ENSG'])
	
	# Table with all transcriptomes
	X = pd.DataFrame(data={'ENSG' : sorted(E)})
	
	for (_, XX) in Progress()(sorted(P2XX.items())) :
		for (n, data) in enumerate(XX) :
			X = pd.merge(X, data, on='ENSG', how='outer')
	
	return { 'FI' : FI, 'X' : X }


def main() :
	
	bundle = get_combined_table()
	
	bundle['script'] = THIS
	
	print("Writing files...")
	
	# Write into a pickle file
	pickle.dump(
		bundle,
		open(OFILE['combined'].format(ext="pkl"), 'wb')
	)
	
	# Write into a txt file
	if PARAM['dump txt'] :
		bundle['X'].to_csv(OFILE['combined'].format(ext="txt"), sep='\t', index=False)


## ==================== ENTRY :

if (__name__ == "__main__") :
	main()
