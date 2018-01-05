import random
import pickle
import math
import networkx as nx
import matplotlib.pyplot as plt
from scipy.constants import golden as phi
from itertools import chain

IFILE = {
	'GO=>Info' : "../OUTPUT/0_go-graph/UV/go-graph.pkl",
}


class gograph :

	#def edge_type(G, e) :
		#return set(G.get_edge_data(*e).keys())

	#for e in GO.edges() :
		#print(edge_type(GO, e))

	def __init__(self, graph_file) :
		self.G = pickle.load(open(graph_file, 'rb'))

		R = {
			# http://amigo.geneontology.org/amigo/term/GO:0008150
			'biological_process',
			
			# http://amigo.geneontology.org/amigo/term/GO:0005575
			'cellular_component',
			
			# http://amigo.geneontology.org/amigo/term/GO:0003674
			'molecular_function',
		}
		#
		R = {
			n : [i for i in self.G.nodes() if (self.G.nodes[i]['name'] == n)]
			for n in R
		}
		#
		assert(all((len(I) == 1) for (n, I) in R.items()))
		#
		self.R = { n : I[0] for (n, I) in R.items() }
		#
		#print("Roots:", self.R)

	def subgraph(self, i, depth=8) :
		I = [i]
		g = nx.MultiDiGraph()
		while (I and (depth > 0)) :
			E = [(i, j) for i in I for j in sorted(self.G.predecessors(i))]
			E = [(i, j) for (i, j) in E if ("obsolete" not in self.G.nodes[j]['name'])]
			I = [j for (i, j) in E if (j not in g.nodes())]
			g.add_edges_from(E)
			depth -= 1

		return nx.MultiDiGraph(nx.subgraph(self.G, g.nodes()))


GO = gograph(IFILE['GO=>Info'])

root = 'GO:0006281' # DNA repair
go = GO.subgraph(root)
print(nx.info(go))



#def subtree(G, i, g, depth=-1) :
	#if (depth == 0) : return
	#J = [j for j in sorted(G.predecessors(i)) if (j not in g.nodes)]
	#for j in J : g.add_edge(i, j)
	#for j in J : subtree(G, j, g, depth-1)
#for depth in range(1, 15) :
	#go = nx.DiGraph()
	#subtree(GO, R['biological_process'], go, depth=depth)
	#print("Number of nodes:", go.number_of_nodes())

def abbr(t) :
	D = { 
		"negative regulation" : "-regu",
		"positive regulation" : "+regu",
		"replication" : "repl",
		"regulation" : "regu",
		"involved in" : "in",
		"synthesis" : "syn",
		"double" : "dbl",
		"single" : "sgl",
		"error" : "err",
	}
	
	for (S, s) in sorted(D.items(), key=(lambda x : -len(x[0]))) :
		t = t.replace(S, s)
	
	return t

from networkx.drawing.nx_agraph import graphviz_layout

## write dot file to use with graphviz
## run "dot -Tpng test.dot >test.png"

pos = graphviz_layout(go, prog='dot', root=root)
pos = { i : (x, -y) for (i, (x, y)) in pos.items() }
plt.figure(figsize=(5, 15), dpi=150)
plt.axis('off')
nx.draw_networkx_nodes(go, pos, node_size=1)
nx.draw_networkx_edges(go, pos, width=0.01, arrows=False, edge_color='y')
#nx.draw_networkx_labels(go, pos, font_size=2)
for (n, (x, y)) in pos.items() :
	al = { 'horizontalalignment' : 'left', 'verticalalignment' : 'bottom' }
	plt.text(x, y, abbr(go.nodes[n]['name'])[0:45] + " ({})".format(n[3:]), fontsize=3, rotation=60, **al)
#plt.show()
#plt.title('')
plt.savefig('test_gotree.eps')

