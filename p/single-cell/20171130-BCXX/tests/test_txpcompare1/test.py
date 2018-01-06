
# RA, 2018-01-06

## ================== IMPORTS :


import os
import pickle
import pandas
import matplotlib.pyplot as plt

## ==================== INPUT :

IFILE = {
	'GO=>CI'   : "../../OUTPUT/0_go2ci/UV/go2ci.pkl",
	'TXP'      : "R00006_TNBC_Go-Scoring-Nolim.csv",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

# See below

## ==================== PARAM :

## ====================== AUX :

## ====================== (!) :

## ===================== DATA :


#[ TXP FILE ]#

txp = pandas.read_csv(IFILE['TXP'])
txp_goid = txp['GO-ID'].tolist()
txp_nsvp = txp['NSV Percentage'].tolist()
txp_ngen = txp['Gene Numbers'].tolist()

#[ MY FILE ]#

# Clustering indices data bundle
CI_data = pickle.load(open(IFILE['GO=>CI'], 'rb'))

# GO2E : GO ID --> Clustering index
GO2CI = CI_data['GO2CI']

# GO2E : GO ID --> [ENSG IDs]
GO2E = CI_data['GO2E']

# GO2T : GO ID --> GO category name
GO2T = CI_data['GO2T']

# N2CI : size of GO term --> [clustering indices]
N2CI = CI_data['N2CI']

# GO2WQ : GO ID --> windowed quantile
GO2WQ = CI_data['GO2WQ']


## ===================== WORK :

plt.scatter(txp_nsvp, [GO2CI[go] for go in txp_goid])
plt.xlabel("Negative silhouette value fraction -- TXP")
plt.ylabel("Average sign of silhouette value -- RA")
plt.savefig("ci.png")
plt.show()

plt.scatter(txp_ngen, [len(GO2E[go]) for go in txp_goid])
plt.xlabel("GO category size -- TXP")
plt.ylabel("GO category size -- RA")
plt.savefig("go.png")
plt.show()

