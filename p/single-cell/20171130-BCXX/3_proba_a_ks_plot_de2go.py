
# RA, 2017-12-06

## ================== IMPORTS :

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

## ==================== INPUT :

input_file_e2ks = "OUTPUT/3_proba_a/e2ks.pkl"

# ENSG ---> GO IDs map
input_file_e2go = "OUTPUT/0_e2go/e2go.txt"

## =================== OUTPUT :

output_file_deplot = "OUTPUT/3_proba_a/de2go_{zoom}.{extension}"
output_file_gorank = "OUTPUT/3_proba_a/gorank.{extension}"

## =================== PARAMS :

num_biomart_parallel_queries = 10
num_biomart_ids_per_query = 100

## ===================== WORK :

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
	print(input_file_e2go, "not found. Try running 0_e2go.py.")
	exit()

# Sorted (differential-expression, gene ID) pairs 
DE2E = sorted((de, e) for (e, de) in E2DE.items())
# Just the DE in the same order
DE = [de for (de, _) in DE2E]

# All known GO IDs
GO = list(set(chain.from_iterable(E2GO.values())))
#
print("GO count:", len(GO))
# 
# Map GO ID --> Count
GO2C = { go : np.zeros(len(DE2E)) for go in GO }
#
for (n, (de, e)) in enumerate(DE2E) :
	for go in E2GO.get(e, []) :
		GO2C[go][n] += 1

GO_C = [ (go, np.cumsum(C)) for (go, C) in GO2C.items() ]
del GO2C

# Sort by decreasing number of occurrence
GO_C = sorted(GO_C, key=(lambda go_C : -max(go_C[1])))

#

plt.clf()
plt.plot([max(C) for (go, C) in GO_C])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("GO ID rank")
plt.ylabel("Number of occurrences in BC data")
plt.savefig(output_file_gorank.format(extension="png"))
plt.savefig(output_file_gorank.format(extension="eps"))

#

#print([(go, x[0:10]) for (go, x) in list(GO2C.items())[0:2]])

# Pick the GO IDs with the most frequent occurrence
GO_C = GO_C[:400]

# Normalize
GO_C = [ (go, C / max(C)) for (go, C) in GO_C ]

GO_C = sorted(GO_C, key=(lambda go_C : -np.mean(go_C[1])))
GO_C = GO_C[:30]

#

plt.clf()

for (_, C) in GO_C :
	plt.plot(DE, C, '-')

plt.xlabel("Differential expression cut-off")
plt.ylabel("Fraction of genes involved")
#plt.xscale("log")
#plt.yscale("log")

# Mechanism + number of associated genes
#L = [(m + " ({})".format(len(M2E[m]))) for m in M]
plt.legend([go for (go, C) in GO_C], prop={'size': 6})

plt.savefig(output_file_deplot.format(zoom="full", extension="png"))
plt.savefig(output_file_deplot.format(zoom="full", extension="eps"))
