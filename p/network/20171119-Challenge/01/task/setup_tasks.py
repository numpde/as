
import testgraph as testgraph

def write_boundary_matrices(G, verbose=False, output_path="./") :
	from itertools import combinations as subcliques
	import networkx as nx
	import gc as gc
	import os as os
	import re as re
	
	def DIAGNOSTIC(*params) :
		if verbose :
			print(*params)

	DIAGNOSTIC("Converting maximal cliques to sorted tuples")
	C = [tuple(sorted(mc)) for mc in nx.find_cliques(G)]
	
	# Each maximal clique is represented in C as a sorted tuple
	
	DIAGNOSTIC("Number of maximal cliques: {} ({}M)".format(len(C), round(len(C) / 1e6)))
	
	# Path where boundary operator matrices will be stored
	if not os.path.exists(output_path) : os.makedirs(output_path)
	
	# Remove existing files
	for f in os.listdir(output_path) :
		if re.search("[0-9]+.txt", f) :
			os.remove(os.path.join(output_path, f))
	
	# Collect the filenames of the matrices
	filenames = []
	
	# (k-1)-chain group
	Sl = dict()
	
	# Iterate over the dimension
	for k in range(0, max(len(s) for s in C)) :
		
		DIAGNOSTIC("Computing the {}-chain group".format(k))
		
		# Get all (k+1)-cliques, i.e. k-simplices, from max cliques mc
		Sk = set(c for mc in C for c in subcliques(mc, k+1))
		# Check that each simplex is in increasing order
		assert(all((list(s) == sorted(s)) for s in Sk))
		# Assign an ID to each simplex, in lexicographic order
		# (This ordering makes subsequent computations faster)
		Sk = dict(zip(sorted(Sk), range(0, len(Sk))))
		
		# Sk is now a representation of the k-chain group
		
		DIAGNOSTIC("{}-chain group rank: {}".format(k, len(Sk)))
		
		# Write boundary matrix to file
		filename = os.path.join(output_path, "{}.txt".format(k))
		filenames.append(filename)
		
		with open(filename, "w") as f :

			DIAGNOSTIC("Writing matrix to", filename)
			
			for ks in Sk.keys() :
				print(' '.join(str(Sl[s]) for s in subcliques(ks, k) if s), file=f)
			
			f.flush()

		Sl = Sk
		gc.collect()
	
	return filenames



def main() :
	import subprocess
	import sys
	
	output_path = sys.argv[1]
	worker_file = sys.argv[2]
	
	assert(output_path)
	assert(worker_file)
	
	for f in write_boundary_matrices(testgraph.make(), output_path=output_path) :
		res = subprocess.run([worker_file, "", "/dev/null"], input=f, 
					   stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
		assert(res.returncode == 0)

		rank = int(res.stdout.splitlines()[0])
			
		print(f, rank)


if (__name__ == "__main__") :
    main()
