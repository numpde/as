
# RA, 2018-01-10

## ================== IMPORTS :

import os
import pickle
import pandas
import numpy as np
import matplotlib.pyplot as plt

## ==================== INPUT :

IFILE = {
	'GO=>CI'       : "../../OUTPUT/0_go2ci/UV/go2ci.pkl",
	'GO graph'     : "../../OUTPUT/0_go-graph/UV/go-graph.pkl",

	'TXP GO names' : "../../ORIGINALS/txp/GOID_to_Name_Mapping.csv",
}

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

## =================== OUTPUT :

OFILE = {
	'diff-names' : "./diff-names.txt",
}

## ===================== DATA :



#[ Compare GO names ]#

txp = pandas.read_csv(IFILE['TXP GO names'])
GO2T_txp = { go : name.lstrip() for (go, name) in zip(txp['GOID'], txp['GOName']) }

go_graph = pickle.load(open(IFILE['GO graph'], 'rb'))
GO2T = { go : go_graph.nodes[go]['name'] for go in go_graph.nodes() }

GO = (set(GO2T.keys()) | set(GO2T_txp.keys()))

with open(OFILE['diff-names'], 'w') as f :

	print("GO ID", "OBO", "TXP", sep='\t', file=f)
	for go in sorted(GO) :
		a = (GO2T    ).get(go, None)
		b = (GO2T_txp).get(go, None)
		
		if not (a == b) :
			print(go, a, b, sep='\t', file=f)
