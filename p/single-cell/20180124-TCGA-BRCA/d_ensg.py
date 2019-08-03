
# RA, 2018-01-25


## ================== IMPORTS :

import sys
import pickle
import os.path
import pandas as pd

from time        import sleep
from itertools   import chain
from joblib      import Parallel, delayed
from progressbar import ProgressBar as Progress

# https://pypi.python.org/pypi/biomart/0.8.0
from biomart     import BiomartServer


## ==================== INPUT :

IFILE = {
	# Extract the list of relevant genes from here
	'TCGA' : "OUTPUT/c_make_table/UV/tcga-brca-fpkm.pkl",
}


## =================== OUTPUT :

OFILE = {
	# Save the results here
	'dump' : "OUTPUT/d_ensg/ensg-info.txt",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)


## =================== PARAMS :

# Biomart root URL
biomart_url = "http://asia.ensembl.org/biomart"

# Query parameters
num_biomart_parallel_queries = 10
num_biomart_ids_per_query = 100
num_biomart_max_retry = 5

TESTMODE = ("TEST" in sys.argv)

if TESTMODE:
	num_biomart_ids_per_query = 1

## ===================== WORK :


# Fetch 'attributes' keyed by the first entry of 'attributes'
def biomart_this(dataset, attributes, E, retry=0) :
	
	try :
		
		response = dataset.search(
			{
				'filters' : { attributes[0] : E }, 
				'attributes': attributes 
			}, 
			header=0
		).content.decode("utf-8")
		
		return [r.split('\t') for r in response.splitlines()]
	
	except Exception as e :
		
		# print("Error on try {}: {}".format(retry, e))
		
		if (retry >= num_biomart_max_retry) : raise e
		
		sleep(retry)
		
		return biomart_this(dataset, attributes, E, retry + 1)


# Download ENSG info from ensembl into a dataframe
def download(E, attributes) :

	# Get biomart record by ENSG ID
	# Q&A: https://www.biostars.org/p/3570/

	# See also
	# https://www.ncbi.nlm.nih.gov/guide/howto/find-func-gene/

	# "The Ensembl genome database project"
	# Nucleic Acids Res. 2002 Jan 1; 30(1): 38–41.
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC99161/

	# Drop the version extension, make unique
	E = set(e[0:15] for e in E)

	# Take only ENSG IDs
	E = sorted(e for e in E if e.startswith("ENSG0"))

	if TESTMODE : 
		E = E[0:5]

	print("Downloading info for {} ENSGs...".format(len(E)))

	# Connect to biomart
	biomart = BiomartServer(biomart_url)
	sapiens = biomart.datasets['hsapiens_gene_ensembl']

	#sapiens.show_filters()
	#sapiens.show_attributes()

	# Partition the list L into chunks of size sz
	# https://stackoverflow.com/a/312466/3609568
	def partition(L, sz):
		L = list(L)
		return [ L[x:(x+sz)] for x in range(0, len(L), sz) ]

	# Download ENSG info chunkwise
	# (biomart is queried via an URL that can't be too long)
	data = list(chain.from_iterable(
		Parallel(n_jobs=num_biomart_parallel_queries, batch_size=1)(
			delayed(biomart_this)(sapiens, attributes, E_part) 
			for E_part in Progress()(partition(E, num_biomart_ids_per_query))
		)
	))

	#dump_data_file = "d_ensg_received_data.log"
	#print("Dumping response to", dump_data_file)
	#with open(dump_data_file, 'w') as f :
		#f.write('\n'.join(data))

	# Make sure all data records have the expected number of entries
	#for row in data: row.extend([None] * (len(attributes) - len(row)))
	assert(all((len(row) == len(attributes)) for row in data))

	print("Creating data frame...")

	# Convert to pandas dataframe, with 'attributes' as column names
	df = pd.DataFrame(data=dict(zip(attributes, list(zip(*data))))).sort_values(attributes[0])
	del data

	print("Merging by ENSG...")

	# Make list of strings unique, filter nontrivial, sort, concatenate
	def cat(L) : return '|'.join(sorted(set(x for x in L if x)))
	# Group by ENSG
	df = df.groupby(by=attributes[0], as_index=False).agg(cat)
	# Rearrange columns to desired order
	df = df[attributes]
	
	return df


def main() :
	
	if os.path.isfile(OFILE['dump']) :
		print("Note:", OFILE['dump'], "already exists. Exiting.")
		return
	
	# Dataframe
	df = download(
		# List of relevant ENSG IDs
		pickle.load(open(IFILE['TCGA'], "rb"))['X']['ENSG'],
		
		# The first entry must be 'ensembl_gene_id'
		[ 'ensembl_gene_id', 'hgnc_symbol', 'go_id' ],
	)

	# Save to text file
	df.to_csv(OFILE['dump'], sep='\t', index=False)


if (__name__ == "__main__") :
	main()
