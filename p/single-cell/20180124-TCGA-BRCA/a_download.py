	
# RA, 2018-01-24

# Run as
#    python3 a*.py DOWNLOAD > tmp.log
#    python3 a*.py SORTCASE


## ================== IMPORTS :

import os
import sys
import json
import hashlib
import urllib.request

from pandas import read_table as pandas_read
from time import sleep
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

from shutil import copyfile
from glob import glob
from itertools import chain
from collections import defaultdict


## ==================== INPUT :

IFILE = {
	# Before download
	'manifest'   : "ORIGINALS/TCGA-BRCA-01/gdc_manifest.2018-01-24T03_39_35.725692.txt",
	'files-json' : "ORIGINALS/TCGA-BRCA-01/files.2018-01-24T04_54_05.741627.json",
}


## =================== OUTPUT :

OFILE = {
	# Download files here
	'files' : "ORIGINALS/TCGA-BRCA-01/UV/{uuid}/{name}",
	
	# Deposit a copy here in the sorting phase
	'casefile' : "OUTPUT/a_cases/UV/{caseid}{sep}{name}",
}


## ==================== PARAM :

PARAM = {
	# Number of parallel HTTP requests
	'requests' : 12,
	
	# https://docs.gdc.cancer.gov/API/Users_Guide/Python_Examples/
	'url-data' : "https://api.gdc.cancer.gov/data/{uuid}",
	'url-info' : "https://api.gdc.cancer.gov/files/{uuid}",
	
	# Persistence upon error
	'max retry' : 10,
}


## ====================== AUX :

# https://stackoverflow.com/questions/107705/disable-output-buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 1)


## ===================== WORK :

def download(uuid, name, md5, retry=0) :
	
	if retry : sleep(retry)
	
	try :
		
		fout = OFILE['files'].format(uuid=uuid, name=name)
		
		def md5_ok() :
			with open(fout, 'rb') as f :
				return (md5 == hashlib.md5(f.read()).hexdigest())
		
		do_download = not (os.path.isfile(fout) and md5_ok())
		
		if do_download :
			
			print("Downloading (attempt {}): {}".format(retry, uuid))
			
			url = PARAM['url-data'].format(uuid=uuid)
			
			with urllib.request.urlopen(url) as response :
				data = response.read()
				
			os.makedirs(os.path.dirname(fout), exist_ok=True)
			
			with open(fout, 'wb') as f :
				f.write(data)
		
		return (uuid, retry, md5_ok())
		
	except Exception as e :
		
		print("Error (attempt {}): {}".format(retry, e))
		
		if (retry >= PARAM['max retry']) :
			raise e
		
		return download(uuid, name, md5, retry + 1)


def DOWNLOAD() :

	mani = pandas_read(IFILE['manifest'], index_col=False, sep='\t')

	OK = Parallel(n_jobs=PARAM['requests'], batch_size=5, verbose=0)(
		delayed(download)(*q)
		for q in Progress()(sorted(zip(mani['id'], mani['filename'], mani['md5'])))
	)

	print("md5 success: {}/{} files.".format(sum(ok for (_, _, ok) in OK), len(mani.index)))


def SORTCASE() :
	
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
	
	files = glob(OFILE['files'].format(uuid="*", name="*"))
	
	for (c, F) in Progress()(C.items()) :
		# Loop over case files
		for name in F :
			f_src = [f for f in files if f.endswith(name)]
			if not f_src :
				print("File not found:", name)
				continue
			
			assert(len(f_src) == 1)
			f_src = f_src[0]
			
			f_dst = OFILE['casefile'].format(caseid=c, sep='/', name=name)
			os.makedirs(os.path.dirname(f_dst), exist_ok=True)
			
			copyfile(f_src, f_dst)
	
	copyfile(IFILE['files-json'], OFILE['casefile'].format(caseid='', sep='', name='meta.json'))


## ==================== ENTRY :

if (__name__ == "__main__") :
	if ("DOWNLOAD" in sys.argv) : DOWNLOAD()
	if ("SORTCASE" in sys.argv) : SORTCASE()

