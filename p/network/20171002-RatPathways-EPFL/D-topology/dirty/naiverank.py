
# RA, 2017-11-06

# import importlib.util as iu
# spec = iu.spec_from_file_location("topology", "../topology.py")
# topology = iu.module_from_spec(spec)
# spec.loader.exec_module(topology)

import networkx as nx
import numpy as np
import itertools
from scipy.sparse import lil_matrix

G = nx.fast_gnp_random_graph(10, 0.5, seed=0)

# All max cliques
C = [tuple(sorted(c)) for c in nx.find_cliques(G)]

# S[k] will hold all k-simplices
# S[k][s] is the ID of simplex s
S = []
for k in range(0, max(len(s) for s in C)) :
	# Get all (k+1)-cliques, i.e. k-simplices, from max cliques mc
	Sk = set(c for mc in C for c in itertools.combinations(mc, k + 1))
	# Check that each simplex is in increasing order
	assert (all((list(s) == sorted(s)) for s in Sk))
	# Assign an ID to each simplex, in lexicographic order
	S.append(dict(zip(sorted(Sk), range(0, len(Sk)))))
	
	print("Number of {}-simplices: {}".format(k, len(Sk)))

# D[k] is the boundary operator
#      from the k complex
#      to the k-1 complex
D = [None for _ in S]

# D[0] maps to zero by definition
D[0] = lil_matrix((1, len(S[0])))

# Parents
P = [dict((sk, []) for sk in Sk) for Sk in S]

# Construct D[1], D[2], ...
for k in range(1, len(S)) :
	D[k] = lil_matrix((len(S[k - 1]), len(S[k])))
	SIGN = np.asmatrix([(-1) ** i for i in range(0, k + 1)]).transpose()
	
	for (ks, j) in S[k].items() :
		for s in itertools.combinations(ks, k) :
			P[k - 1][s].append(ks)
		
		# Indices of all (k-1)-subsimplices s of the k-simplex ks
		I = [S[k - 1][s] for s in sorted(itertools.combinations(ks, k))]
		D[k][I, j] = SIGN
		
for (k, d) in enumerate(D) :
	print("D[{}] has shape {}".format(k, d.shape))

# Check that D[k-1] * D[k] is zero
assert(all((0 == np.dot(D[k-1], D[k]).count_nonzero()) for k in range(1, len(D))))


# Compute rank using matrix SVD
rk = [np.linalg.matrix_rank(d.todense()) for d in D] + [0]
# Compute dimker using rank-nullity theorem
ns = [(d.shape[1] - rk[n]) for (n, d) in enumerate(D)] + [0]

print("rk:", rk)
print("ns:", ns)


H = [None for _ in S]
H[0] = nx.Graph(G)
for k in range(1, len(S)) :
	H[k] = nx.Graph()
	for (e, L) in P[k - 1].items() :
		H[k].add_edges_from(itertools.combinations(L, 2))

#for h in H : print(h.edges())

for h in H :
	print(list(h.degree()))
	T = nx.minimum_spanning_tree(h)
	print(len(T.edges()), h.number_of_nodes())