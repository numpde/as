
# RA, 2018-01-29

# Run as
#    python3 f*.py


## ================== IMPORTS :

import os
import pickle
import inspect
import pandas as pd

## ==================== INPUT :

IFILE = {
	# "Hallmarks" gene sets
	# http://software.broadinstitute.org/gsea/msigdb/collections.jsp
	'hallmarks' : "ORIGINALS/GeneSets/broadinstitute.org/V/h.all.v6.1.symbols.gmt",
	
	# GO annotations
	'go->symb' : "ORIGINALS/GO/UV/go2symb.txt",
	'go->name' : "ORIGINALS/GO/UV/go2name.txt",
}


## =================== OUTPUT :

OFILE = {
	# Selected gene subsets
	'subsets' : "OUTPUT/f_gene_subsets/subsets.pkl",
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)

## ==================== PARAM :

PARAM = {
	# GO terms of interest
	'GO' : {
		"GO:0001525", # angiogenesis
		"GO:0006281", # DNA-repair
		#"GO:0006955", # immune response
		"GO:0007049", # cell cycle
		"GO:0016477", # cell migration
		#"GO:0004984", # olfactory receptor activity
		#"GO:0045596", # negative regulation of cell differentiation
		#"GO:0045597", # positive regulation of cell differentiation
		"GO:0000723", # telomere maintenance
		
		"GO:0007173", # epidermal growth factor receptor signaling pathway
		"GO:0035004", # phosphatidylinositol 3-kinase activity
		"GO:0051019", # mitogen-activated protein kinase binding
		
		# Estrogen-related GO terms
		"GO:0030284", # estrogen receptor activity
		#"GO:0030520", # intracellular estrogen receptor signaling pathway
		#"GO:0030331", # estrogen receptor binding
		#"GO:0038049", # transcription factor activity, ligand-activated RNA polymerase II transcription factor binding
		
		# HER2
		"GO:0038128", # ERBB2 signaling pathway (ERBB2 = ENSG00000141736)
	},
}

## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Log which files are opened
def logged_open(filename, mode='r', *argv, **kwargs) :
	print("({}):\t{}".format(mode, filename))
	return open(filename, mode, *argv, **kwargs)


## ====================== (!) :

def geneset(set_id, set_info, set_genes) :
	return { 'id' : set_id, 'info' : set_info, 'set' : set(set_genes) }


## ===================== WORK :

def hallmarks() :

	with logged_open(IFILE['hallmarks'], mode='r') as f :
		S = [
			geneset(l[0], l[1], l[2:])
			for l in [l.rstrip().split('\t') for l in f.readlines()] 
		]
	
	return S


def gocategories() :
	GO = pd.read_table(IFILE['go->symb'], index_col=0).join(pd.read_table(IFILE['go->name'], index_col=0), how='outer')
	
	S = [
		geneset(go, GO['GO name'][go], GO['Primary symbol'][go].split('|'))
		for go in PARAM['GO']
	]
		
	return S


def main() :
	
	S = gocategories() + hallmarks()
	
	pickle.dump(
		{
			'S'      : S,
			'script' : THIS,
			'param'  : PARAM,
		},
		logged_open(OFILE['subsets'], 'wb')
	)


## ==================== ENTRY :

if (__name__ == "__main__") :
	main()

