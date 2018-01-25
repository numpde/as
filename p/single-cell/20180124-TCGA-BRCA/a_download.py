	
# RA, 2018-01-24

# Run as
#    python3 a_download.py > tmp.log

## ================== IMPORTS :

import os
import sys
import hashlib
import urllib.request

from pandas import read_table as pandas_read
from time import sleep
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

## ==================== INPUT :

IFILE = {
	'manifest' : "ORIGINALS/TCGA-BRCA-01/gdc_manifest.2018-01-24T03_39_35.725692.txt",
}

## =================== OUTPUT :

OFILE = {
	'data' : "ORIGINALS/TCGA-BRCA-01/UV/{uuid}/{name}",
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


def job(uuid, name, md5, retry=0) :
	
	if retry : sleep(retry)
	
	try :
		
		fout = OFILE['data'].format(uuid=uuid, name=name)
		
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
		
		return job(uuid, name, md5, retry+1)


mani = pandas_read(IFILE['manifest'], index_col=False, sep='\t')

OK = Parallel(n_jobs=PARAM['requests'], batch_size=5, verbose=0)(
	delayed(job)(*q)
	for q in Progress()(sorted(zip(mani['id'], mani['filename'], mani['md5'])))
)

print("md5 success: {}/{} files.".format(sum(ok for (_, _, ok) in OK), len(mani.index)))
