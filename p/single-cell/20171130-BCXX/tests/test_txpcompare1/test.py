
# RA, 2018-01-06

## ================== IMPORTS :

import os
import pickle
import pandas
import numpy as np
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

## ===================== DATA :


#[ TXP FILE ]#

txp = pandas.read_csv(IFILE['TXP'])
txp_go = txp['GO-ID'].tolist()
txp_ci = txp['NSV Percentage'].tolist()
txp_ng = txp['Gene Numbers'].tolist()

#[ MY FILE ]#

# Clustering indices data bundle
CI_data = pickle.load(open(IFILE['GO=>CI'], 'rb'))

# GO2E : GO ID --> Clustering index
GO2CI = CI_data['GO2CI']

# GO2E : GO ID --> [ENSG IDs]
GO2E = CI_data['GO2E']

# GO2T : GO ID --> GO category name
GO2T = CI_data['GO2T']


## ===================== WORK :

txp_ci = np.asarray(txp_ci)
ra_ci  = np.asarray([GO2CI[go] for go in txp_go])

txp_ng = np.asarray(txp_ng)
ra_ng  = np.asarray([len(GO2E[go]) for go in txp_go])

A = np.asarray([(abs(x - y) <= 5) for (x, y) in zip(txp_ng, ra_ng)])
B = np.logical_not(A)

h = [
	plt.scatter(txp_ci[B], ra_ci[B], color='r'),
	plt.scatter(txp_ci[A], ra_ci[A], color='b'),
]
h.append(
	plt.plot(plt.xlim(), 1-2*np.asarray(plt.xlim()), '--k')[0]
)
plt.xlabel("Negative silhouette value fraction -- TXP")
plt.ylabel("Average sign of silhouette value -- RA")
plt.legend(h, ["size diff > 5", "size diff <= 5", "expected"], loc='lower left')
plt.savefig("ci.png")
plt.show()

h = [
	plt.scatter(txp_ng[B], ra_ng[B], color='r'),
	plt.scatter(txp_ng[A], ra_ng[A], color='b'),
]
h.append(
	plt.plot(plt.xlim(), plt.xlim(), '--k')[0]
)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("GO category size -- TXP")
plt.ylabel("GO category size -- RA")
plt.legend(h, ["size diff > 5", "size diff <= 5", "expected"], loc='lower left')
plt.savefig("go.png")
plt.show()

