
import numpy as np

input_file_hgnc = "ORIGINALS/UV/HGNC-Approved-20171205.txt"

input_file_func = { 
	"DNA repair" : "ORIGINALS/byfunction/dna-repair.txt",
	"Differentiation" : "ORIGINALS/byfunction/differentiation.txt",
}

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


for (func, filename) in input_file_func.items() :
	print("Mechanism:", func)
	
	X = np.genfromtxt(filename, delimiter='\t', dtype=None, comments='\r')

	U = [u[len("UniProtKB:"):] for u in utf(X[:, 0])]
	E = [U2E[u] for u in U if (u in U2E)]
	print("Mapped: {} genes. Not mapped: {}".format(len(E), [u for u in U if (u not in U2E)]))

