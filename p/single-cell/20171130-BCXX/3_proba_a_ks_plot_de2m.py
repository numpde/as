
# RA, 2017-12-06

### IMPORTS -- #

import pickle
import numpy as np
import matplotlib.pyplot as plt

from itertools import chain

### INPUT ---- #

input_file_e2ks = "OUTPUT/3_proba_a/UV/gene-ks.pkl"

input_file_hgnc = "ORIGINALS/UV/HGNC-Approved-20171205.txt"

input_file_mech = { 
	"Apoptosis"       : "ORIGINALS/byfunction/GO/apoptosis.txt",
	"DNA repair"      : "ORIGINALS/byfunction/GO/dna-repair.txt",
	"Differentiation" : "ORIGINALS/byfunction/GO/differentiation.txt",
	"BRCA"            : "ORIGINALS/byfunction/GO/brca.txt",
	"Cell migration"  : "ORIGINALS/byfunction/GO/migration.txt",
	"Angiogenesis"    : "ORIGINALS/byfunction/GO/angiogenesis.txt",
	"Cell cycle"      : "ORIGINALS/byfunction/GO/cellcycle.txt",
	"Immune response" : "ORIGINALS/byfunction/GO/immuneresponse.txt"
}

### OUTPUT --- #

output_file_deplot = "OUTPUT/3_proba_a/de2m_{zoom}.{extension}"

### MEAT ----- #

# Get the [ENSG ID] <--> [UniProt ID] mappings

def utf(S) : return [s.decode("utf-8") for s in S]

def map_uniprot_ensg(filename) :
	X = np.genfromtxt(filename, delimiter='\t', dtype=None, names=True, comments='\r')
	u = [n for n in X.dtype.names if ("UniProt" in n)].pop()
	e = [n for n in X.dtype.names if ("Ensembl" in n)].pop()
	UE = [(u, e) for (u, e) in zip(utf(X[u]), utf(X[e])) if (u and e)]
	U2E = dict((u, e) for (u, e) in UE)
	E2U = dict((e, u) for (u, e) in UE)
	return (U2E, E2U)

(U2E, E2U) = map_uniprot_ensg(input_file_hgnc)

# Sanity check
assert(E2U[U2E["Q9H2K8"]] == "Q9H2K8")

# Load the Differential Expression data (DE)
# For a gene e, E2DE[e] is a measure of DE
E2DE = pickle.load(open(input_file_e2ks, "rb"))['E2KS']

# Load the mechanism-to-genes maps
M2E = dict()
for (mech, filename) in input_file_mech.items() :
	print("Mechanism:", mech)
	
	Y = np.genfromtxt(filename, delimiter='\t', dtype=None, comments='\r')

	U = [u[len("UniProtKB:"):] for u in utf(Y[:, 0])]
	E = [U2E[u] for u in U if (u in U2E)]
	print("Mapped: {} genes. Not mapped: {}".format(len(E), [u for u in U if (u not in U2E)]))

	# List of genes as ENSG IDs
	M2E[mech] = [e for e in E if (e in E2DE)]

M2E["Other"] = list(set(E2DE.keys()) - set(chain.from_iterable(M2E.values())))
M2E["All"] = E2DE.keys()

# Sorted (differential-expression, gene ID) pairs 
DE2E = sorted((de[1], e) for (e, de) in E2DE.items())
# Just the DE in the same order
DE = [de for (de, _) in DE2E]

# Mechanism-to-count
M2C = {
	m : np.cumsum([(e in E) for (_, e) in DE2E])
	for (m, E) in M2E.items()
}

# Normalize
M2C = { m : C / max(C) for (m, C) in M2C.items() }

# List of mechanisms sorted by average involvement
M = reversed(sorted((np.mean(C), m) for (m, C) in M2C.items()))
M = [m for (_, m) in M]

for m in M :
	plt.plot(DE, M2C[m], '-')

plt.xlabel("Differential expression cut-off")
plt.ylabel("Fraction of genes involved")

# Mechanism + number of associated genes
L = [(m + " ({})".format(len(M2E[m]))) for m in M]
plt.legend(L)

plt.savefig(output_file_deplot.format(zoom="full", extension="png"))
plt.savefig(output_file_deplot.format(zoom="full", extension="eps"))
