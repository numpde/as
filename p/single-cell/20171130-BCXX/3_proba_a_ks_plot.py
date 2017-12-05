
# RA, 2017-12-05

### IMPORTS -- #

import pickle
import numpy as np
import matplotlib.pyplot as plt

### INPUT ---- #

input_file_e2ks = "OUTPUT/3_proba_a/e2ks.pkl"

input_file_hgnc = "ORIGINALS/UV/HGNC-Approved-20171205.txt"

input_file_func = { 
	"Apoptosis"       : "ORIGINALS/byfunction/apoptosis.txt",
	"DNA repair"      : "ORIGINALS/byfunction/dna-repair.txt",
	"Differentiation" : "ORIGINALS/byfunction/differentiation.txt",
	"BRCA"            : "ORIGINALS/byfunction/brca.txt",
	"Cell migration"  : "ORIGINALS/byfunction/migration.txt",
	"Angiogenesis"    : "ORIGINALS/byfunction/angiogenesis.txt",
	"Cell cycle"      : "ORIGINALS/byfunction/cellcycle.txt",
	"Immune response" : "ORIGINALS/byfunction/immuneresponse.txt"
}

### OUTPUT --- #

output_file_deplot = "OUTPUT/3_proba_a/de_{zoom}.{extension}"

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
E2DE = pickle.load(open(input_file_e2ks, "rb"))['E2DE']

# Cycle through the cell mechanisms

M2DE = dict()
for (mech, filename) in input_file_func.items() :
	print("Mechanism:", mech)
	
	X = np.genfromtxt(filename, delimiter='\t', dtype=None, comments='\r')

	U = [u[len("UniProtKB:"):] for u in utf(X[:, 0])]
	E = [U2E[u] for u in U if (u in U2E)]
	print("Mapped: {} genes. Not mapped: {}".format(len(E), [u for u in U if (u not in U2E)]))

	de = list(sorted(E2DE[e] for e in E if (e in E2DE)))
	M2DE[mech] = de

# Mechanisms
M = list(sorted((np.mean(de), m) for (m, de) in M2DE.items()))
M = [m for (_, m) in M]

for m in M :
	plt.plot(M2DE[m], np.linspace(0.0, 1.0, len(M2DE[m]) + 1)[1:], '-')

plt.xlabel("Differential expression (average interbatch KS distance)")
plt.ylabel("Gene relative rank in mechanism")
#plt.xscale("log")
#plt.yscale("log")

# Mechanism + number of associated genes
L = [(m + " ({})".format(len(M2DE[m]))) for m in M]
plt.legend(L)

plt.savefig(output_file_deplot.format(zoom="full", extension="png"))
plt.savefig(output_file_deplot.format(zoom="full", extension="eps"))

plt.xlim((plt.xlim()[0], np.mean(plt.xlim())))
plt.ylim((plt.ylim()[0], np.mean(plt.ylim())))

plt.savefig(output_file_deplot.format(zoom="zoom", extension="png"))
plt.savefig(output_file_deplot.format(zoom="zoom", extension="eps"))
