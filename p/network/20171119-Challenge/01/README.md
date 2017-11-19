The purpose of this challenge is to improve the rank computation routine for graph-topological boundary operators.

A reference implementation is provided in the subfolder **reference**. It reads a file name from standard input and prints the result to standard output. The file is expected to contain a binary matrix in the format

	137 148 406
	140 149 450
	75 81 225
	72 77 154
	181 194 435
	...

where a number **i** the **j**-th row indicates that the (i,j) component is nonzero. The output of the routine is the rank of the matrix over Z/2Z.

The file **task/testgraph.py** should contain a routine **make()** that returns the graph used for testing, for example

	def make() : 
		import networkx as nx
		# Erdos-Renyi random graph
		n = 40
		p = 0.7
		G = nx.gnp_random_graph(n, p, seed=0)
		return G

From this graph, the matrices of boundary operators of every dimension will be constructed.

Run

	make setup

to compile the reference implementation, create the matrices, write them to disk and run the reference implementation on them. Depending on the size of the testgraph, this may take a while. The list of created files and the output of the reference implementation are contained in the file **tasks.txt**, e.g. for the above graph, 

	./task/UV/0.txt 0
	./task/UV/1.txt 40
	./task/UV/2.txt 550
	./task/UV/3.txt 3438
	./task/UV/4.txt 10842
	./task/UV/5.txt 18037
	./task/UV/6.txt 15845
	./task/UV/7.txt 6670
	./task/UV/8.txt 1078

Name the executable of your implementation **rank** and put it into the folder **your**. Then run

	make check

to compare the results. The output should look as follows, where "expected" is the result of the reference implementation and "received" is the result of your implementation:

	Testfile: ./task/UV/0.txt
	Expected: 0
	Received: 0

	Testfile: ./task/UV/1.txt
	Expected: 40
	Received: 40

	Testfile: ./task/UV/2.txt
	Expected: 550
	Received: 550

	Testfile: ./task/UV/3.txt
	Expected: 3438
	Received: 3438

	Testfile: ./task/UV/4.txt
	Expected: 10842
	Received: 10842

	Testfile: ./task/UV/5.txt
	Expected: 18037
	Received: 18037

	Testfile: ./task/UV/6.txt
	Expected: 15845
	Received: 15845

	Testfile: ./task/UV/7.txt
	Expected: 6670
	Received: 6670

	Testfile: ./task/UV/8.txt
	Expected: 1078
	Received: 1078
