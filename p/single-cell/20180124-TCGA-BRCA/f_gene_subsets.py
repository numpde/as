
# RA, 2018-01-29

# Run as
#    python3 f*.py


## ================== IMPORTS :

import os
import pickle
import inspect
import pandas as pd
from progressbar import ProgressBar as Progress

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
		
		# From B_go_by_wq, gene numbers are approximate
		"GO:0005132", # type I interferon receptor binding (6)
		#"GO:0015874", # norepinephrine transport (3)
		# "GO:0005333", # norepinephrine transmembrane transporter activity (3)
		"GO:0008188", # neuropeptide receptor activity (12)
		"GO:0004956", # prostaglandin D receptor activity (2)
		"GO:0014050", # negative regulation of glutamate secretion (5)
		"GO:0033038", # bitter taste receptor activity (19)
		"GO:1904315", # transmitter-gated ion channel activity involved in regulation of postsynaptic membrane potential (3)
		"GO:0042571", # immunoglobulin complex, circulating (30)
		"GO:0086006", # voltage-gated sodium channel activity involved in cardiac muscle cell action potential (5)
		"GO:1901374", # acetate ester transport (2)
		#"GO:0005277", # acetylcholine transmembrane transporter activity (2)
		"GO:1901078", # negative regulation of relaxation of muscle (3)
		"GO:0060078", # regulation of postsynaptic membrane potential (45)
		"GO:0002323", # natural killer cell activation involved in immune response (8)
		"GO:1905144", # response to acetylcholine (3)
		"GO:0035377", # transepithelial water transport (2)
		"GO:0007196", # adenylate cyclase-inhibiting G-protein coupled glutamate receptor signaling pathway (8)
		"GO:0005344", # oxygen carrier activity (12)
		"GO:0060371", # regulation of atrial cardiac muscle cell membrane depolarization (6)
		"GO:0030718", # germ-line stem cell population maintenance (3)
		"GO:0008527", # taste receptor activity (14)
		"GO:0015747", # urate transport (5)
		"GO:0070698", # type I activin receptor binding (2)
		"GO:0042693", # muscle cell fate commitment (4)
		"GO:0072014", # proximal tubule development (3)
		"GO:0086010", # membrane depolarization during action potential (26)
		"GO:0050907", # detection of chemical stimulus involved in sensory perception (28)
		"GO:0005249", # voltage-gated potassium channel activity (67)
		"GO:0001518", # voltage-gated sodium channel complex (13)
		"GO:0005261", # cation channel activity (29)
		"GO:0030280", # structural constituent of epidermis (10)
		"GO:0005245", # voltage-gated calcium channel activity (36)
		"GO:0019373", # epoxygenase P450 pathway (17)
		"GO:0008395", # steroid hydroxylase activity (20)
		"GO:0006855", # drug transmembrane transport (22)
		"GO:0007218", # neuropeptide signaling pathway (80)
		"GO:0002377", # immunoglobulin production (55)
		"GO:0004970", # ionotropic glutamate receptor activity (17)
		"GO:2000330", # positive regulation of T-helper 17 cell lineage commitment (4)
		"GO:0015671", # oxygen transport (13)
		"GO:0017144", # drug metabolic process (23)
		"GO:0006828", # manganese ion transport (10)
		"GO:0007608", # sensory perception of smell (174)
	},
	
	# https://www.biostars.org/p/77590/
	# Changed ORC6L to ORC6 (https://en.wikipedia.org/wiki/ORC6, 2019-08-02)
	'PAM50' : "ACTR3B, ANLN, BAG1, BCL2, BIRC5, BLVRA, CCNB1, CCNE1, CDC20, CDC6, CDH3, CENPF, CEP55, CXXC5, EGFR, ERBB2, ESR1, EXO1, FGFR4, FOXA1, FOXC1, GPR160, GRB7, KIF2C, KRT14, KRT17, KRT5, MAPT, MDM2, MELK, MIA, MKI67, MLPH, MMP11, MYBL2, MYC, NAT1, NDC80, NUF2, ORC6, PGR, PHGDH, PTTG1, RRM2, SFRP1, SLC39A6, TMEM45B, TYMS, UBE2C, UBE2T",
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

	for go in PARAM['GO']:
		if go not in GO.index:
			print("{} unknown".format(go))
	
	S = [
		geneset(go, GO['GO name'][go], GO['Primary symbol'][go].split('|'))
		for go in PARAM['GO']
	]
		
	return S


def pam50() :
	s = geneset("PAM50", "PAM50", PARAM['PAM50'].replace(", ", "\t").split())
	
	return [ s ]

def main() :
	
	S = pam50() + gocategories() + hallmarks()
	
	for s in S :
		print("{} ({})".format(s['id'], len(s['set'])))
	
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

