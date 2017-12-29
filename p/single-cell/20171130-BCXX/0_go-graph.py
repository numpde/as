
# RA, 2017-12-25

import obonet
import networkx
import pickle

pickle.dump(
	obonet.read_obo(
		"http://purl.obolibrary.org/obo/go.obo"
	),
	open(
		"OUTPUT/0_go-graph/UV/go-graph.pkl",
		'wb'
	)
)
