	
# RA, 2018-01-24

# Run as
#    python3 b1*.py


## ================== IMPORTS :

import os
import sys
import math
import pickle
import inspect
import xmltodict
import pandas as pd

from glob import glob

from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

from itertools import chain
from collections import OrderedDict


## ==================== INPUT :

IFILE = {
	# Pattern for finding casefiles
	'casefile' : "OUTPUT/a_cases/UV/{caseid}/{name}",
}


## =================== OUTPUT :

OFILE = {
	'clinical' : "OUTPUT/b1_clinical/tcga-clinical.{ext}",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Number of parallel computing processes
	'#proc' :  int(TESTMODE) or min(12, math.ceil(cpu_count() / 1.5)),
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Change key to conform with the workflow
def rekey(k) :
	if k.endswith("followup_barcode") : return "followup_barcode"
	if k.endswith("patient_barcode") : return "patient_barcode"
	if k.endswith("patient_uuid") : return "patient_uuid"
	return k


## ===================== WORK :

def clinical(x) :
	if (type(x) is str) :
		# Interpret x as filename
		with open(x, mode='r') as f :
			return clinical(xmltodict.parse(f.read())['brca:tcga_bcr']['brca:patient'])
	
	if (type(x) is list) :
		return list(chain.from_iterable(clinical(y) for y in x))
	
	# Expect x to be a parsed xml node
	if not (type(x) is OrderedDict) :
		print(type(x))
		print(x)
	
	assert(type(x) is OrderedDict)
	
	records = []
	
	record = { }
	#
	for (k, v) in x.items() :
		try : 
			record[rekey(v['@preferred_name'] or k)] = v['#text']
		except :
			pass
	
	records.append(pd.Series(data=record))
	
	for (k, y) in x.items() :
		if k.endswith("follow_ups") :
			try :
				for z in y.values() :
					records.extend(clinical(z))
			except :
				pass
	
	return records


def CLINICAL() :
	
	print("Parsing files...")
	
	# List of clinical files
	files = glob(IFILE['casefile'].format(caseid="*", name="*clinical*.xml"))
	
	if TESTMODE : files = files[0:10]
	
	C = list(chain.from_iterable(
		Parallel(n_jobs=PARAM['#proc'], batch_size=5)(
			delayed(clinical)(f)
			for f in Progress()(files)
		)
	))

	
	print("Merging records...")
	
	DATA = {
		# Patient info (temporary container)
		'patient_barcode' : [],
		
		# Follow-up info (temporary container)
		'followup_barcode' : [],
	}
	
	for c in Progress()(C) :
		for k in DATA.keys() :
			if (k in c.index) :
				DATA[k].append(c.rename(c[k]))
				c = None
				break
		assert(not c), "Record has to be either a 'patient' or a 'follow-up'"
	
	for (k, data) in DATA.items() :
		DATA[k] = pd.DataFrame().join(data, how='outer')
	
	
	print("Saving...")
	
	# Save data as text files
	for (k, data) in DATA.items() :
		data.to_csv(OFILE['clinical'].format(ext=(k + ".txt")), sep='\t')
	
	# Save data in binary
	DATA['script'] = THIS
	DATA['param'] = PARAM
	pickle.dump(
		DATA,
		open(OFILE['clinical'].format(ext="pkl"), 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	CLINICAL()

