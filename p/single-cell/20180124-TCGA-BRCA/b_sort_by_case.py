
# RA, 2018-01-24

# Run as
#    python3 b*.py

## ================== IMPORTS :

import gzip
import sys
import os
import json

from shutil import copyfile
from glob import glob
from itertools import chain
from collections import defaultdict

## ==================== INPUT :

IFILE = {
	'files-json' : "ORIGINALS/TCGA-BRCA-01/files.2018-01-24T04_54_05.741627.json",
	'files'      : glob("ORIGINALS/TCGA-BRCA-01/UV/*/*"),
}

## =================== OUTPUT :

OFILE = {
	'casefile' : "OUTPUT/b_cases/UV/{caseid}{sep}{name}",
}


## ==================== PARAM :

PARAM = {
}


## ===================== WORK :


def sort_by_case() :
	
	L = json.load(open(IFILE['files-json']))
	assert(all((len(rec['cases']) == 1) for rec in L))
	
	C = defaultdict(list)
	for rec in L :
		C[rec['cases'][0]['case_id']].append( rec['file_name'] )
	C = dict(C)
	
	# Test
	#for rec in L :
		#if ('0a2a3529-f645-4967-9a58-89ee20b8bb62' == rec['cases'][0]['case_id'] ) :
			#print(rec)
			#exit()
	
	for (c, F) in C.items() :
		# Loop over case files
		for name in F :
			f_src = [f for f in IFILE['files'] if f.endswith(name)]
			if not f_src :
				print("File not found:", name)
				continue
			
			assert(len(f_src) == 1)
			f_src = f_src[0]
			
			f_dst = OFILE['casefile'].format(caseid=c, sep='/', name=name)
			os.makedirs(os.path.dirname(f_dst), exist_ok=True)
			
			copyfile(f_src, f_dst)
	
	copyfile(IFILE['files-json'], OFILE['casefile'].format(caseid='', sep='', name='meta.json'))


## ===================== MAIN :

if (__name__ == "__main__") :
	sort_by_case()
