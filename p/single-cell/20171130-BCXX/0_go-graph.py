
# RA, 2017-12-25

import obonet
import networkx
import pickle

go_obo_url = "http://purl.obolibrary.org/obo/go.obo"
graph_file = "OUTPUT/0_go-graph/UV/go-graph.pkl"

pickle.dump(
	obonet.read_obo(
		go_obo_url
	),
	open(
		graph_file,
		'wb'
	)
)
