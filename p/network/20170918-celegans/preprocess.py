#!/usr/bin/python3

# RA, 2017-09-18

# Part I
# Read the /C. elegans/ connectome from the excel file
# 	NeuronConnect.xls
# downloaded from
# 	http://www.wormatlas.org/neuronalwiring.html
# 	http://www.wormatlas.org/images/NeuronConnect.xls
#
# Part II
# Convert it to graph objects and save to disk
#
# Part III
# Test-read from disk


###########################
print('\n', "Part I", '\n')
###########################

filename_in = "./NeuronConnect.xls"
url = "http://www.wormatlas.org/images/NeuronConnect.xls"


import os.path
if not os.path.isfile(filename_in) :
	# Download the file
	print("File", filename_in, "not found. Downloading.")
	import urllib.request
	urllib.request.urlretrieve(url, filename_in)
	print("Done.")
	print(" ")

assert(os.path.isfile(filename_in))


from xlrd import open_workbook
# https://stackoverflow.com/questions/22169325/read-excel-file-in-python

wb = open_workbook('NeuronConnect.xls')

synapses = []
for sheet in wb.sheets() :
	for nrow in range(1, sheet.nrows):
		synapse = [x.value for x in sheet.row(nrow)]

		# Neuromuscular junction
		is_nmj = (synapse[2] == 'NMJ')

		synapses.append(synapse)

neurons = []
neurons += [s[0] for s in synapses]
neurons += [s[1] for s in synapses]
neurons = set(neurons)

# Now have the set of unique /neurons/ :
print("Neurons:", neurons)
# and the list of /synapses/ :
print("First synapse:", synapses[0])
print("Last synapse: ", synapses[-1])


############################
print('\n', "Part Ia", '\n')
############################

# Do some preprocessing and consistency checks

# Uppercase all
synapses = [[s[0].upper(), s[1].upper(), s[2], s[3]] for s in synapses]

# Check that every EJ synapse is present symmetrically but not repeated
S_EJ = dict([])
for s in synapses :
	if (s[2] != 'EJ') : continue
	s1 = [s[1], s[0], s[2], s[3]]
	assert(s1 in synapses)

	k = (s[0], s[1])
	v = s[3]
	assert(not k in S_EJ)
	S_EJ[k] = v

# Check that NMJ synapses are present only once
S_NMJ = dict([])
for s in synapses :
	if (s[2] != 'NMJ') : continue
	assert(s[0] != 'NMJ')
	assert(s[1] == 'NMJ')

	k = (s[0], s[1])
	v = s[3]
	S_NMJ[k] = v

# Collapse multiply defined synapses
S_send = dict([]) # Will be used for the graph
S_recv = dict([]) # Used only for checking
for s in synapses :
	c = s[2]
	if (c == 'EJ') : continue
	if (c == 'NMJ') : continue
	assert(c in ['R', 'S', 'Rp', 'Sp'])

	if (c in ['S', 'Sp']) :
		k = (s[0], s[1])
		if (not k in S_send) : S_send[k] = 0
		S_send[k] += s[3]

	if (c in ['R', 'Rp']) :
		k = (s[1], s[0])
		if (not k in S_recv) : S_recv[k] = 0
		S_recv[k] += s[3]

assert(len(set(S_recv.keys()) - set(S_send.keys())) == 0)
assert(len(set(S_send.keys()) - set(S_recv.keys())) == 0)

for (k, v) in S_send.items() :
	assert(S_recv[k] == v)

del S_recv
del synapses

############################
print('\n', "Part II", '\n')
############################

# CONSTRUCT GRAPHS

import networkx as nx

# Empty directed graph for chemical synapses
G = nx.DiGraph()
# Empty directed graph for EJ synapses (all double!)
E = nx.DiGraph()
# Empty directed graph for NMJ synapses
M = nx.DiGraph()

# Populate the graph with nodes
G.add_nodes_from(neurons)
E.add_nodes_from(neurons)
M.add_nodes_from(neurons)


print("Graph nodes:", G.nodes())

for (k, v) in S_send.items() :
	G.add_edge(k[0], k[1], weight=v)

for (k, v) in S_EJ.items() :
	E.add_edge(k[0], k[1], weight=v)

for (k, v) in S_NMJ.items() :
	M.add_edge(k[0], k[1], weight=v)

# SAVE GRAPHS IN PICKLE FORMAT

import pickle

filename_pkl = "./celegans.pkl"

with open(filename_pkl, "wb") as f :
	pickle.dump(G, f, protocol=pickle.HIGHEST_PROTOCOL)
	pickle.dump(E, f, protocol=pickle.HIGHEST_PROTOCOL)
	pickle.dump(M, f, protocol=pickle.HIGHEST_PROTOCOL)

print("Saved to file:", filename_pkl)


# SAVE THE G GRAPH IN TXT FORMAT

filename_txt = "./celegans_{}.txt"

first_label = 0
G0 = nx.convert_node_labels_to_integers(G, first_label=first_label)

with open(filename_txt.format("head"), "w") as f :
	def P(X) : print(X, file=f)
	P("# Directed *multigraph* of chemical synapses of C. Elegans")
	P("# From: http://www.wormatlas.org/images/NeuronConnect.xls")
	P("# First node: {}".format(first_label))
	P("# number of nodes = {}".format(G0.number_of_nodes()))
	W = [data['weight'] for (a, b, data) in G.edges(data=True)]
	for w in W : assert(w.is_integer() and (w >= 0))
	P("# number of edges = {}".format(int(sum(W))))

with open(filename_txt.format("data"), "w") as f :
	def P(X) : print(X, file=f)
	for (a, b, data) in G0.edges(data=True) :
		for _ in range(int(data['weight'])) :
			P("{} {}".format(a, b))

# PLOT

import matplotlib.pyplot as plt
plt.ion()

plt.figure()
nx.draw(G)
plt.show()

plt.figure()
nx.draw(E)
plt.show()

plt.figure()
nx.draw(M)
plt.show()


#############################
print('\n', "Part III", '\n')
#############################

# Test read pickle file


del G, E, M

with open(filename_pkl, "rb") as f :
	G = pickle.load(f)
	E = pickle.load(f)
	M = pickle.load(f)

print("Nodes:", len(G.nodes()))
print("S/R edges:", len(G.edges()))
print("EJ edges:", len(E.edges()))
print("NMJ edges:", len(M.edges()))

# Count directed synapses with multiplicity
print("Number of directed synapses:", sum(data['weight'] for (a, b, data) in G.edges(data=True)))

input("Press enter to finish...")


