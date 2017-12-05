
# RA, 2017-12-02

# (In progress)

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
}

### OUTPUT --- #

output_file_deplot = "OUTPUT/3_proba_a/de.{extension}"

### MEAT ----- #

# Get the ENSG <--> UniProt mappings

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

# 

E2KS = pickle.load(open(input_file_e2ks, "rb"))['E2KS']

# Cycle through the cell mechanisms

M2KS = dict()
for (mech, filename) in input_file_func.items() :
	print("Mechanism:", mech)
	
	X = np.genfromtxt(filename, delimiter='\t', dtype=None, comments='\r')

	U = [u[len("UniProtKB:"):] for u in utf(X[:, 0])]
	E = [U2E[u] for u in U if (u in U2E)]
	print("Mapped: {} genes. Not mapped: {}".format(len(E), [u for u in U if (u not in U2E)]))

	ks = list(reversed(sorted(E2KS[e] for e in E if (e in E2KS))))
	M2KS[mech] = ks


M = list(sorted((np.mean(ks), k) for (k, ks) in M2KS.items()))
print(M)
M = [k for (_, k) in M]

for m in M :
	plt.plot(range(1, 1 + len(M2KS[m])), M2KS[m], '--.')

plt.xlabel("Involved gene")
plt.ylabel("Differential expression")
plt.xscale("log")

plt.savefig(output_file_deplot.format(extension="png") + "-noleg.png")

plt.legend(M)

plt.savefig(output_file_deplot.format(extension="png"))

plt.show()
