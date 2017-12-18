
# RA, 2017-12-06

## ================== IMPORTS :

import pickle
import os.path

from itertools   import chain
from progressbar import ProgressBar as Progress
from joblib      import Parallel, delayed

#https://pypi.python.org/pypi/biomart/0.8.0
from biomart     import BiomartServer

## ==================== INPUT :

IFILE = {
	# Extract the list of relevant genes from here
	'BC selected' : "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl",
}

## =================== OUTPUT :

OFILE = {
	# ENSG ID ---> GO ID associations
	'E2GO' : "OUTPUT/0_e2go/e2go.txt",
	'All GO' : "OUTPUT/0_e2go/go.txt",
}


## =================== PARAMS :

# 
biomart_url = "http://grch37.ensembl.org/biomart"

num_biomart_parallel_queries = 10
num_biomart_ids_per_query = 100

## ===================== WORK :

# [ PART 1 ]

if os.path.isfile(OFILE['E2GO']) :
	
	print(OFILE['E2GO'], "already exists.")
	
	# Read the ENSG-GO associations from file
	with open(OFILE['E2GO'], 'r') as f :
		E_GO = [L.rstrip().split('\t') for L in f]
		E2GO = { e_go[0] : e_go[1:] for e_go in E_GO }
		del E_GO
		
	# E2GO[e] is now a list of GO IDs associated to ENSG ID e
	
	print("E2GO:", list(E2GO.items())[0:10])

else :
	# Download ENSG-GO associations
	
	# Get biomart record by ENSG ID
	# Q&A: https://www.biostars.org/p/3570/
	
	# See also
	# https://www.ncbi.nlm.nih.gov/guide/howto/find-func-gene/
	
	# "The Ensembl genome database project"
	# Nucleic Acids Res. 2002 Jan 1; 30(1): 38–41.
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC99161/
	
	# List of relevant ENSG IDs
	E = pickle.load(open(IFILE['BC selected'], "rb"))['gene_id']
	
	print("Downloading ENSG-GO associations from biomart...")
	
	# Connect to biomart
	biomart = BiomartServer(biomart_url)
	sapiens = biomart.datasets['hsapiens_gene_ensembl']

	#sapiens.show_filters()
	#sapiens.show_attributes()
		
	def biomart_e2go(E) :
		
		response = sapiens.search(
			{
				'filters' : { 'ensembl_gene_id' : E },
				'attributes': [ 'ensembl_gene_id', 'go_id' ]
			}, 
			header=0
		).content.decode("utf-8")
		
		return [r.split() for r in response.splitlines()]

	# Partition the list L into chunks of size sz
	# https://stackoverflow.com/a/312466/3609568
	def partition(L, sz):
		L = list(L)
		return [L[x:(x+sz)] for x in range(0, len(L), sz)]

	# Download ENSG-GO associations from biomart chunkwise
	# (biomart is queried via an URL that can't be too long)
	E_GO = list(chain.from_iterable(
		Parallel(n_jobs = num_biomart_parallel_queries)(
			delayed(biomart_e2go)(E_part) 
			for E_part in Progress()(partition(E, num_biomart_ids_per_query))
		)
	))
	
	assert(len(E_GO)), "E_GO is empty"
	
	# Filter the ENSG that have any associated GO
	E_GO = [e_go for e_go in E_GO if (len(e_go) == 2)]
	
	E2GO = dict()
	for (e, go) in E_GO : E2GO.setdefault(e, []).append(go)
	
	with open(OFILE['E2GO'], 'w') as f :
		for (e, GO) in E2GO.items() :
			print('\t'.join([e] + GO), file=f)

# [ PART 2 ]

with open(OFILE['All GO'], 'w') as f :
	for go in sorted(set(chain.from_iterable(E2GO.values()))) :
		print(go, file=f)
