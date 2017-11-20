
import time
class Timer(object):
	def __init__(self, name=None):
		self.name = name

	def __enter__(self):
		self.tstart = time.time()

	def __exit__(self, type, value, traceback):
		print('[{}], elapsed: {}'.format(self.name, (time.time() - self.tstart)))


def betti_bin_cpp(C, verbose=False) :
	from itertools import combinations as subcliques
	import subprocess as subprocess
	import tempfile as tempfile
	import gc as gc
	import os as os
	
	def DIAGNOSTIC(*params) :
		if verbose :
			print(*params)

	if ((type(C) is list) and all((mc == tuple(sorted(mc))) for mc in C)) :
		DIAGNOSTIC("Using the provided maximal cliques list")
		pass
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
		
		if (k == 0) :
			Sl = Sk
			gc.collect()
			continue
		
		# CALL THE EXTERNAL ROUTINE
		
		with tempfile.NamedTemporaryFile(mode="w") as f :
		
			filename = f.name
			
			for ks in Sk.keys() :
				print(' '.join(str(Sl[s]) for s in subcliques(ks, k) if s), file=f)
			
			f.seek(0)
			f.flush()
			
			#with Timer("subprocess") :
			res = subprocess.run(["./a.out", "", "/dev/null"], input=filename, stdout=subprocess.PIPE, universal_newlines=True)
			assert(res.returncode == 0)
			rank = int(res.stdout.splitlines()[0])
		
		DIAGNOSTIC("rank = ", rank)
		
		B[k]   -= rank
		B[k-1] -= rank
		
		if (not verbose) : 
			S = None
			Sl = Sk
			gc.collect()
			continue
		
		# CROSS-CHECK
		
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
		
		rank_ref = 0
		
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
			
			rank_ref += 1
		
		DIAGNOSTIC("rank_ref = ", rank_ref)
		
		assert(rank == rank_ref)
		
		gc.collect()
	
	while B and (B[-1] == 0) : B.pop()
	
	return B

def test() :
	import networkx as nx
	from test_topology_localcopy import betti_bin as betti_bin_ref
	
	G = nx.gnp_random_graph(40, 0.77, seed=0)

	with Timer('c++') :
		print(betti_bin_cpp(nx.find_cliques(G)))
	
	#with Timer('ref') :
		#print(betti_bin_ref(nx.find_cliques(G)))
	

test()
