
# Compute the Betti numbers of a graph
# (from its maximal cliques)
#
# C is an unused networkx.find_cliques(G) of some graph
#
# RA, 2017-11-03 (CC-BY-4.0)
#
# Ref:
#   A. Zomorodian, Computational topology (Notes), 2009
#   http://www.ams.org/meetings/short-courses/zomorodian-notes.pdf
#
# Gist:
#   https://gist.github.com/numpde/16f3a22e352dc43dc01614b50b74645b
#
def betti(C, verbose = False) :
	
	import itertools
	import numpy as np
	import networkx as nx
	from scipy.sparse import lil_matrix
	
	def DIAGNOSTIC(*params) :
		if verbose :
			print(*params)
	
	#
	# 1. Prepare maximal cliques
	#
	
	# Sort each maximal clique, make sure it's a tuple
	# Also, commit C to memory if necessary
	C = [tuple(sorted(c)) for c in C]
	
	DIAGNOSTIC("Number of maximal cliques: {} ({}M)".format(len(C), round(len(C) / 1e6)))
	
	#
	# 2. Enumerate all simplices
	#
	
	# S[k] will hold all k-simplices
	# S[k][s] is the ID of simplex s
	S = []
	for k in range(0, max(len(s) for s in C)) :
		# Get all (k+1)-cliques, i.e. k-simplices, from max cliques mc
		Sk = set(c for mc in C for c in itertools.combinations(mc, k+1))
		# Check that each simplex is in increasing order
		assert(all((list(s) == sorted(s)) for s in Sk))
		# Assign an ID to each simplex, in lexicographic order
		S.append(dict(zip(sorted(Sk), range(0, len(Sk)))))

		DIAGNOSTIC("Number of {}-simplices: {}".format(k, len(Sk)))
	
	# The cliques are redundant now
	del C
		
	# Euler characteristic
	ec = sum(((-1)**k * len(S[k])) for k in range(0, len(S)))
	
	DIAGNOSTIC("Euler characteristic:", ec)
		
	#
	# 3. Construct the boundary operators
	#
	
	# D[k] is the boundary operator
	#	  from the k complex
	#	  to the k-1 complex
	D = [None for _ in S]
		
	# D[0] maps to zero by definition
	D[0] = lil_matrix( (1, len(S[0])) )

	# Construct D[1], D[2], ...
	for k in range(1, len(S)) :
		D[k] = lil_matrix( (len(S[k-1]), len(S[k])) )
		SIGN = np.asmatrix([(-1)**i for i in range(0, k+1)]).transpose()

		for (ks, j) in S[k].items() :
			# Indices of all (k-1)-subsimplices s of the k-simplex ks
			I = [S[k-1][s] for s in sorted(itertools.combinations(ks, k))]
			D[k][I, j] = SIGN

	# The simplices are redundant now
	del S
	
	for (k, d) in enumerate(D) :
		DIAGNOSTIC("D[{}] has shape {}".format(k, d.shape))
	
	# Check that D[k-1] * D[k] is zero
	assert(all((0 == np.dot(D[k-1], D[k]).count_nonzero()) for k in range(1, len(D))))
		
	#
	# 4. Compute rank and dimker of the boundary operators
	#
	
	# Compute rank using matrix SVD
	rk = [np.linalg.matrix_rank(d.todense()) for d in D] + [0]
	# Compute dimker using rank-nullity theorem
	ns = [(d.shape[1] - rk[n]) for (n, d) in enumerate(D)] + [0]

	# The boundary operators are redundant now
	del D

	DIAGNOSTIC("rk:", rk)
	DIAGNOSTIC("ns:", ns)

	#
	# 5. Infer the Betti numbers
	#
	
	# Betti numbers
	# B[0] is the number of connected components
	B = [(n - r) for (n, r) in zip(ns[:-1], rk[1:])]

	# Remove trailing zeros
	while B and (B[-1] == 0) : B.pop()

	DIAGNOSTIC("Betti numbers:", B)
	
	# Check: Euler-Poincare formula (see Eqn 16 in [Zomorodian])
	epc = sum(((-1)**k * B[k]) for k in range(0, len(B)))
	if (ec != epc) : print("Warning: ec = {} and epc = {} should be equal".format(ec, epc))

	return B


# Compute the Betti numbers of a graph over Z/2Z
# 
# If G is a netwokx graph then
# C = networkx.find_cliques(G)
# (enumerates maximal cliques)
#
# RA, 2017-11-08 (CC-BY-4.0)
#
# Ref: 
#   A. Zomorodian, Computational topology (Notes), 2009
#   http://www.ams.org/meetings/short-courses/zomorodian-notes.pdf
#
# The computation of the rank is adapted from
#   https://triangleinequality.wordpress.com/2014/01/23/computing-homology/
# see also
#   https://gist.github.com/numpde/9584779ad235c6ee19be7a6bb87e8af5
#
# Gist for this function:
#   https://gist.github.com/numpde/a3f3476e601041da3cb968f31e95409d
#
def betti_bin(C, verbose=False) :
	from itertools import combinations as subcliques
	
	def DIAGNOSTIC(*params) :
		if verbose :
			print(*params)
	
	if all((mc == tuple(sorted(mc))) for mc in C) :
		DIAGNOSTIC("Using the provided maximal cliques list")
	else :
		DIAGNOSTIC("Converting maximal cliques to sorted tuples")
		C = [tuple(sorted(mc)) for mc in C]
		
	# Each maximal clique is represented in C as a sorted tuple
	
	DIAGNOSTIC("Number of maximal cliques: {} ({}M)".format(len(C), round(len(C) / 1e6)))
	
	# Betti numbers
	B = []
	
	# (k-1)-chain group
	Sl = dict()
	
	# Iterate over the dimension
	for k in range(0, max(len(s) for s in C)) :
		
		# Get all (k+1)-cliques, i.e. k-simplices, from max cliques mc
		Sk = set(c for mc in C for c in subcliques(mc, k+1))
		# Check that each simplex is in increasing order
		assert(all((list(s) == sorted(s)) for s in Sk))
		# Assign an ID to each simplex, in lexicographic order
		# (This ordering makes subsequent computations faster)
		Sk = dict(zip(sorted(Sk), range(0, len(Sk))))
		
		# Sk is now a representation of the k-chain group
		
		DIAGNOSTIC("{}-chain group rank: {}".format(k, len(Sk)))
		
		# "Default" Betti number (if rank is 0)
		B.append(len(Sk))
		
		# J2I is a mapped representation of the boundary operator
		# (Understood to have boolean entries instead of +1 / -1)
		
		J2I = {
			j : set(Sl[s] for s in subcliques(ks, k) if s)
			for (ks, j) in Sk.items()
		}
		
		Sl = Sk
		
		DIAGNOSTIC("Assembled J2I of size {}".format(len(J2I)))
		
		# Transpose J2I
		
		I2J = {
			i : set()
			for I in J2I.values() for i in I
		}
		
		for (j, I) in J2I.items() :
			for i in I :
				I2J[i].add(j)
		
		DIAGNOSTIC("Assembled I2J of size {}".format(len(I2J)))
				
		# Compute the rank of the boundary operator
		# The rank = The # of completed while loops
		# It's usually the most time-consuming part
		while J2I and I2J :
			# Default pivot option
			I = next(iter(J2I.values()))
			J = I2J[next(iter(I))]
			
			# An alternative option
			JA = next(iter(I2J.values()))
			IA = J2I[next(iter(JA))]
			
			# Choose the option with fewer indices
			# (reducing the number of loops below)
			if ((len(IA) + len(JA)) < (len(I) + len(J))) : (I, J) = (IA, JA)

			for i in I : 
				I2J[i] = J.symmetric_difference(I2J[i])
				if not I2J[i] : del I2J[i]

			for j in J : 
				J2I[j] = I.symmetric_difference(J2I[j])
				if not J2I[j] : del J2I[j]
			
			# Let
			#   D_k : C_k --> C_{k-1}
			# be the boundary operator. Recall that
			#   H_k = ker D_k / im D_{k+1}

			# Rank of  D_k  increases, therefore:
			B[k-1] -= 1
			# Nullspace dimension of  D_k  decreases, therefore:
			B[k]   -= 1
			
			#DIAGNOSTIC("Potential pivots left:", len(J2I) + len(I2J))
	
	# Remove trailing zeros
	while B and (B[-1] == 0) : B.pop()
	
	return B
	