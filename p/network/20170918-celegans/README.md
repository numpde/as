Construct the graph of the C. elegans neural connectome.

The graph is constructed from the file

	NeuronConnect.xls

that was obtained (2017-09-18) from 

	http://www.wormatlas.org/neuronalwiring.html
	http://www.wormatlas.org/images/NeuronConnect.xls

Run

	python3 xls2graph.py

to extract the graph from the xls file. 

The xls file lists the synapses:

	Neuron 1	Neuron 2	Type	Nbr
	ADAR	ADAL	EJ	1
	ADFL	ADAL	EJ	1
	ASHL	ADAL	EJ	1
	AVDR	ADAL	EJ	2
	PVQL	ADAL	EJ	1
	ADEL	ADAL	Sp	1
	ADFL	ADAL	Sp	1
	AIAL	ADAL	Sp	1
	AIBL	ADAL	R	1
	AIBR	ADAL	Rp	2
	ASHL	ADAL	Sp	1
	...
	VD12	NMJ	NMJ	13
	VD13	NMJ	NMJ	12

The 3rd column **Type** denotes the type of the synapse,

	S   = Send or output (Neuron 1 pre-synaptic to Neuron 2)
	Sp  = Send-poly (Neuron 1 is pre-synaptic to more than one postsynaptic partner)
	R   = Receive or input (Neuron 1 is post-synaptic to Neuron 2)
	Rp  = Receive-poly (Neuron 1 is one of several post-synaptic partners of Neuron 2)
	EJ  = Electric junction
	NMJ = Neuromuscular junction (only reconstructed NMJ's are represented)

and the 4th column **Nbr** is the number of synapses between the given neuron pair.

The routine **xls2graph.py** extracts the three graphs **G**, **E** and **M**, that contain, respectively, the chemical synapses (of type **S**/**Sp**/**R**/**Rp**), the **EJ** synapses, and the **NMJ** synapses. The graphs are stored in the file

	celegans.pkl

The graphs are weighted directed graphs, where the weight is the multiplicity of the synapse.
Note that the xls file contains all synapses twice, expect those of type **NMJ**.
The graph **E** of electric junctions is a directed graph where the edge **(a, b)** is present whenever **(b, a)** is.

The graph of the chemical synapses is also written to the two files

**celegans_head.txt**:

	# Directed *multigraph* of chemical synapses of C. Elegans
	# From: http://www.wormatlas.org/images/NeuronConnect.xls
	# First node: 0
	# number of nodes = 283
	# number of edges = 6394

that contains meta-data on the connectivity graph,
and

**celegans_data.txt**:

	0 209
	0 245
	1 172
	1 172
	1 172
	1 130
	1 68
	1 68
	2 161
	2 161
	...

that lists the directed edges (synapses) between the nodes (neurons) according to their multiplicity.
