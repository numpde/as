
# RA, 2018-01-18

# Run as
#    python3 A*.py download process

## ================== IMPORTS :

# For downloading and uncompressing files
import urllib.request, gzip

# For command line arguments
import sys

# For making folders and looking up files
import os

# For reading and handling tables
from pandas import read_table

# For handling graphs
import networkx as nx

# For concatenating lists and such
from itertools import chain

# For grouping pairs by first entry
from itertools import groupby
from operator import itemgetter

## ==================== INPUT :

pass


## ==================== PARAM :

PARAM = {
	# GO annotations
	# http://geneontology.org/gene-associations/readme/goa_human.README
	'GOA URL gaf' : "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
}

## =================== OUTPUT :

OFILE = {
	'GOA' : "ORIGINALS/UV/goa_human.gaf",
	
	'synonyms' : "OUTPUT/synonyms.txt",
	'symb->go' : "OUTPUT/symb2go.txt",
	'go->symb' : "OUTPUT/go2symb.txt",
	'syn2symb' : "OUTPUT/syn2symb.txt",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)


## ===================== WORK :


# Download GO annotation file
def download() :
	
	def get(remote, local) :
		# https://stackoverflow.com/a/7244263/3609568
		# Download the zipped file
		with urllib.request.urlopen(remote) as response :
			# Unzip t
			with gzip.GzipFile(fileobj=response) as uncompressed :
				# Open a local file for writing a bytes object
				with open(local, 'wb') as f :
					# Write the unzipped contents
					f.write(uncompressed.read())
	
	get(PARAM['GOA URL gaf'], OFILE['GOA'])


# Read the GO annotations table
# Perform some consistency checks
# Make bipartite graphs 
#  H of primary symbols <--> synonyms
#  G of primary symbols <--> GO terms
# Returns the triple of (set of primary symbols, G, H)
def read_goa(filename) :
	
	# Use pandas to load the GO annotations table
	# http://www.geneontology.org/page/go-annotation-file-format-20
	names = ['DB', 'ID', 'Symbol', 'Q', 'GO', 'DB Ref', 'Evidence', 'With/From', 'Aspect', 'Name', 'Synonyms', 'Type', 'Taxon', 'Date', 'Assigned by', 'Extension', 'Gene product ID']
	data = read_table(filename, index_col=False, comment='!', sep='\t', header=None, nrows=None, names=names)
	
	# Check that each ID has at exactly one (primary) symbol
	# (but some symbols have multiple IDs)
	g = nx.Graph()
	g.add_edges_from(('ID:' + i, 'Sy:' + s) for (i, s) in zip(data['ID'], data['Symbol']))
	assert(all((d == 1) for (n, d) in g.degree() if n.startswith('ID:')))
	del g
	
	# Obtain the symbol synonyms
	SYN = [
		(s, set(syn.split('|')) - set(s))
		for (s, syn) in zip(data['Symbol'], data['Synonyms'])
	]
	
	# Check that each symbol has a unique list of synonyms
	assert(all(
		(1 == len(set(tuple(sorted(s)) for (_, s) in S)))
		for (_, S) in groupby(SYN, itemgetter(0))
	))
	
	# Now can collaps the list to a dict
	SYN = dict(SYN)
	
	# Remove the (primary) symbol from the synonyms
	SYN = {
		s0 : [s for s in syn if (s not in SYN.keys())]
		for (s0, syn) in SYN.items()
	}
	
	# Primary symbols
	S0 = set(SYN.keys())
	# All synonyms (excluding the primary symbols)
	S1 = set(chain.from_iterable(SYN.values()))
	# They should be disjoint by now:
	assert(not (S0 & S1))
	
	# Make a bipartite graph "ID" <--> "Synonyms"
	H = nx.Graph()
	H.add_nodes_from(S0, is_symbol=1, is_synonym=0)
	H.add_nodes_from(S1, is_symbol=0, is_synonym=1)
	H.add_edges_from((i, s) for (i, syn) in SYN.items() for s in syn)
	# There should be no self-loops
	assert(all((a != b) for (a, b) in H.edges()))
	# Each synonym belongs to at least one primary symbol
	for s in S1 : assert(H.degree(s) >= 1)
	
	# Make a bipartite graph "Primary gene symbols" <--> "GO terms"
	# https://networkx.github.io/documentation/stable/reference/algorithms/bipartite.html
	G = nx.Graph()
	assert(S0 == set(data['Symbol']))
	G.add_nodes_from(data['Symbol'], is_symbol=1, is_go=0)
	G.add_nodes_from(data['GO'],     is_symbol=0, is_go=1)
	G.add_edges_from(zip(data['Symbol'], data['GO']))
	
	return (S0, G, H)


def process() :
	
	(S0, G, H) = read_goa(OFILE['GOA'])
	
	# Symbols and synonyms
	
	with open(OFILE['synonyms'], 'w') as f :
		print("Primary symbol", "Non-ambiguous synonyms", "Ambiguous synonyms", sep='\t', file=f)
		for s in sorted(S0) : 
			# Ambiguous synonyms
			a = sorted(n for n in H.neighbors(s) if (H.degree(n) >= 2))
			# Non-ambiguous synonyms
			b = sorted(n for n in H.neighbors(s) if (H.degree(n) == 1))
			print(s, '|'.join(b), '|'.join(a), sep='\t', file=f)

	with open(OFILE['syn2symb'], 'w') as f :
		print("Synonym", "Primary symbols", "Ambiguity", sep='\t', file=f)
		for s in sorted(H.nodes()) :
			if H.nodes[s]['is_synonym'] :
				n = sorted(H.neighbors(s))
				print(s, '|'.join(n), (len(n) - 1), sep='\t', file=f)
	
	# GO terms
	
	with open(OFILE['symb->go'], 'w') as f :
		print("Primary symbol", "GO terms", sep='\t', file=f)
		for s in sorted(S0) : 
			print(s, '|'.join(sorted(G.neighbors(s))), sep='\t', file=f)

	with open(OFILE['go->symb'], 'w') as f :
		print("GO term", "Primary symbol", sep='\t', file=f)
		for s in sorted(set(G.nodes()) - S0) : 
			assert(s.startswith("GO:"))
			print(s, '|'.join(sorted(G.neighbors(s))), sep='\t', file=f)



if (__name__ == "__main__") :

	if ("download" in sys.argv) : download()
	if ("process"  in sys.argv) : process()


