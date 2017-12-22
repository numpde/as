
# RA, 2017-12-02

# Select only the meaningful gene records from the BC data

### IMPORTS -- #

import re
import os
import pickle
import inspect
import numpy as np

### INPUT ---- #

input_file_raw_BC = "ORIGINALS/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
input_file_select = "ORIGINALS/txp/gene-selection.txt"

### OUTPUT --- #

output_file_select = "OUTPUT/0_select/UV/" + os.path.split(input_file_raw_BC)[1] + "-selected.pkl"

### MEAT ----- #

# Load the raw data
X = np.genfromtxt(input_file_raw_BC, names=True, delimiter='\t', dtype=None)

# Feature tags from the raw file 
header = list(X.dtype.names)
#print(header[0:20])

# Filter the BCXX(LN)_XX columns
p = re.compile(r"BC[0-9]+.*_[0-9]+", re.IGNORECASE)
header = [n for n in header if p.match(n)]

# The gene type of each gene
gene_type = [t.decode("utf-8") for t in X['gene_type']]

# The ENSGXXXXXXXXXXX gene id of each gene
gene_id = [t.decode("utf-8")[0:15] for t in X['gene_id']]

# Histogram of gene types
#gene_type_freq = set((t, len([g for g in gene_type if (g == t)])) for t in set(gene_type))
#print("Gene types:", set(gene_type_freq))
## {('misc_RNA', 2034), ('TR_V_gene', 97), ('rRNA', 527), ('IG_V_pseudogene', 187), ('Mt_tRNA', 22), ('sense_overlapping', 202), ('IG_J_pseudogene', 3), ('pseudogene', 13931), ('SPIKE_IN', 3), ('sense_intronic', 742), ('TR_C_gene', 5), ('IG_D_gene', 37), ('TR_V_pseudogene', 27), ('miRNA', 3055), ('TR_J_pseudogene', 4), ('TR_D_gene', 3), ('snRNA', 1916), ('polymorphic_pseudogene', 45), ('IG_V_gene', 138), ('TR_J_gene', 74), ('lincRNA', 7114), ('antisense', 5276), ('IG_C_pseudogene', 9), ('snoRNA', 1457), ('3prime_overlapping_ncrna', 21), ('IG_J_gene', 18), ('protein_coding', 20345), ('IG_C_gene', 14), ('processed_transcript', 515), ('ERCC', 92), ('Mt_rRNA', 2)}

# Make a numpy matrix: (metadata + samples) x genes
X = np.asarray(list(X[n] for n in header))

# Along which axis/dimension of the data are the samples / the genes
(axis_smpl, axis_gene) = (0, 1)

# Consistency check
assert(X.shape[axis_gene] == len(gene_id)), "Gene IDs not in sync"
assert(X.shape[axis_gene] == len(gene_type)), "Gene types not in sync"

#
def delbyid(L, J) :
	return [v for (j, v) in enumerate(L) if (j not in set(J))]

if False :
	# Omit the BC01 cluster
	
	# Cluster membership of each sample (BCXX)
	sample2cluster = [n[0:4] for n in header]
	
	I = [i for (i, c) in enumerate(sample2cluster) if (c == "BC01")]
	X = np.delete(X, I, axis_smpl)
	sample2cluster = [c for (i, c) in enumerate(sample2cluster) if (i not in I)]
	del header # invalidate inconsitent variable

# Omit uninteresting gene records

if False : 
	# Keep only genes of this type:
	T = {"protein_coding"}
	
	# Which genes to omit?
	J = [j for (j, t) in enumerate(gene_type) if (t not in T)]
	X = np.delete(X, J, axis_gene)
	(gene_id, gene_type) = (delbyid(gene_id, J), delbyid(gene_type, J))
	print("Omitted {} genes that are not {}".format(len(J), T))

if True :
	# Use the selection file to filter relevant genes
	Y = np.genfromtxt(input_file_select, names=True, delimiter='\t', dtype=None)
	
	# 15 is the length of the ENSGXXXXXXXXXXX gene_id code
	Y = set([t.decode("utf-8")[0:15] for t in Y['gene_id']])
	
	# Which genes to omit?
	J = [j for (j, i) in enumerate(gene_id) if not (i in Y)]
	X = np.delete(X, J, axis_gene)
	(gene_id, gene_type) = (delbyid(gene_id, J), delbyid(gene_type, J))
	print("Omitted {} genes not listed in {}".format(len(J), input_file_select))

if True :
	# Omit genes that are not expressed
	J = np.nonzero(np.std(X, axis=axis_smpl) == 0)[0].tolist()
	X = np.delete(X, J, axis_gene)
	(gene_id, gene_type) = (delbyid(gene_id, J), delbyid(gene_type, J))
	print("Omitted {} non-expressed genes".format(len(J)))

# Size of the remaining dataset
(n_smpls, n_genes) = (X.shape[axis_smpl], X.shape[axis_gene])
#
print("Got data with {} samples and {} nontrivial genes".format(n_smpls, n_genes))

# Sanity check
assert(n_smpls), "Didn't get any samples"
assert(n_genes), "Didn't get any genes"
assert(X.shape[axis_gene] == len(gene_id)), "Gene IDs not in sync"
assert(X.shape[axis_gene] == len(gene_type)), "Gene types not in sync"

# Group samples by patient BCXX
#
patients = [h[0:4] for h in header]
#
P2SH = {
	p : [(s, h) for (s, h) in enumerate(header) if re.match(p + "_[0-9]+", h)]
	for p in patients 
}

# Groups samples by batch BCXX[LN][_Re]
#
batches = [re.findall("(.*)_[0-9]+", h)[0] for h in header]
#
B2SH = {
	b : [(s, h) for (s, h) in enumerate(header) if re.match(b + "_[0-9]+", h)]
	for b in batches 
}


# To get the gene-wise z-score transform do
#
#	from scipy import stats
#	Z = stats.mstats.zscore(X, axis=axis_smpl)

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
script = inspect.getsource(inspect.getmodule(inspect.currentframe()))

# Save the filtered data
pickle.dump(
	{
		'X'         : X,            # filtered data
		'header'    : header,       # BCXX[LN]_[Re]_XX
		'gene_id'   : gene_id,      # ENSGXXXXXXXXXXX
		'gene_type' : gene_type,    # protein_coding, lincRNA, ...
		'n_smpls'   : n_smpls,      # number of samples
		'n_genes'   : n_genes,      # number of genes
		'axis_smpl' : axis_smpl,    # axis/dimension of samples
		'axis_gene' : axis_gene,    # axis/dimension of genes
		'P2SH'      : P2SH,         # group by patient
		'B2SH'      : B2SH,         # group by batch
		'script'    : script        # this script
	},
	open(output_file_select, "wb")
)
