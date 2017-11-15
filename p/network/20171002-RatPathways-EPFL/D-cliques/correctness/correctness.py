
# RA, 2017-11-15

# Check the computation of Betti numbers
# against an independent implementation

### IMPORTS -- #

import importlib.util as iu
spec = iu.spec_from_file_location("topology", "../topology.py")
topology = iu.module_from_spec(spec)
spec.loader.exec_module(topology)

import networkx as nx
import numpy as np
import os
from subprocess import call
from tempfile import mktemp
from scipy.io import savemat, loadmat
from scipy.sparse import csc_matrix

### MEAT ----- #

### Graph to test implementations

G = nx.fast_gnp_random_graph(30, 0.4)

### Compute the local candidate 

bc = topology.betti(nx.find_cliques(G))

print("Candidate:", bc)

### Invoke MATLAB for reference 

M = csc_matrix(nx.adjacency_matrix(G), dtype='d')

filename = mktemp(".mat")
#filename = "./tmp.mat"
savemat(filename, { 'G' : M })

x = '-nodesktop -nojvm -r "' + "betti('{}'); exit;".format(filename) + '"';
call(["matlab", x + " > /dev/null"])

br = loadmat(filename)['b'].flatten().tolist()
while br and (br[-1] == 0) : br.pop()

os.remove(filename)

print("Reference:", br)
