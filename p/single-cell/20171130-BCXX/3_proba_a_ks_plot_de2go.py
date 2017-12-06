
# RA, 2017-12-06

## ------------------ IMPORTS :

import pickle
import os.path
import itertools
import numpy as np
import matplotlib.pyplot as plt

from itertools   import chain
from progressbar import ProgressBar as Progress
from joblib      import Parallel, delayed

#https://pypi.python.org/pypi/biomart/0.8.0
from biomart     import BiomartServer

## -------------------- INPUT :

input_file_e2ks = "OUTPUT/3_proba_a/e2ks.pkl"

# If the file
input_file_e2go = "OUTPUT/3_proba_a/e2go.txt"
# does not exist, it will be created via biomart
biomart_url = "http://grch37.ensembl.org/biomart"

## ------------------- OUTPUT :

output_file_deplot = "OUTPUT/3_proba_a/de2go_{zoom}.{extension}"

# Will be created if necessary:
output_file_e2go = input_file_e2go

## ------------------- PARAMS :

num_biomart_parallel_queries = 10
num_biomart_ids_per_query = 100

## --------------------- WORK :

# Load the Differential Expression data (DE)
# For a gene e, E2DE[e] is a measure of DE
E2DE = pickle.load(open(input_file_e2ks, "rb"))['E2DE']

# Get ENSG --> GOs associations

if os.path.isfile(input_file_e2go) :
	# Read ENSG-GO associations from file
	
	# E2GO[e] is a list of GO IDs associated to ENSG ID e
	with open(input_file_e2go, 'r') as f :
		E_GO = [L.rstrip().split('\t') for L in f]
		E2GO = { e_go[0] : e_go[1:] for e_go in E_GO }
		del E_GO

else :
	# Download ENSG-GO associations
	
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
	E_GO = list(itertools.chain.from_iterable(
		Parallel(n_jobs = num_biomart_parallel_queries)(
			delayed(biomart_e2go)(E_part) 
			for E_part in Progress()(partition(E2DE.keys(), num_biomart_ids_per_query))
		)
	))
	
	assert(len(E_GO)), "E_GO is empty"
	
	# Filter the ENSG that have any associated GO
	E_GO = [e_go for e_go in E_GO if (len(e_go) == 2)]
	
	E2GO = dict()
	for (e, go) in E_GO : E2GO.setdefault(e, []).append(go)
	
	with open(output_file_e2go, 'w') as f :
		for (e, GO) in E2GO.items() :
			print('\t'.join([e] + GO), file=f)
	
	del E_GO
	del sapiens
	del biomart

# Now E2GO[e] is a list of GO IDs associated to ENSG ID e


# Sorted (differential-expression, gene ID) pairs 
DE2E = sorted((de, e) for (e, de) in E2DE.items())
# Just the DE in the same order
DE = [de for (de, _) in DE2E]

# All known GO IDs
GO = list(set(chain.from_iterable(E2GO.values())))
# 
# Map GO ID --> Count
GO2C = { go : np.zeros(len(DE2E)) for go in GO }
#
for (n, (de, e)) in enumerate(DE2E) :
	for go in E2GO.get(e, []) :
		GO2C[go][n] += 1
#
GO2C = { go : np.cumsum(C) for (go, C) in GO2C.items() }

#print([(go, x[0:10]) for (go, x) in list(GO2C.items())[0:2]])

GO = reversed(sorted((max(C), go) for (go, C) in GO2C.items()))
GO = [go for (_, go) in GO]
GO = GO[0:20]

for go in GO:
	plt.plot(DE, GO2C[go], '-')

#plt.ylabel("Proportion of genes involved")
#plt.xlabel("Differential expression cut-off")
#plt.xscale("log")
#plt.yscale("log")

# Mechanism + number of associated genes
#L = [(m + " ({})".format(len(M2E[m]))) for m in M]
plt.legend(GO, prop={'size': 6})

#plt.savefig(output_file_deplot.format(zoom="full", extension="png"))
#plt.savefig(output_file_deplot.format(zoom="full", extension="eps"))

plt.show()
