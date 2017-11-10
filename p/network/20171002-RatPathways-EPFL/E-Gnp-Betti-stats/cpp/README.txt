The c++ routine rank.cpp computes the Z/2Z rank of a binary matrix. 

It reads a filename from stdin.

The file is expected to be in the format:

	1 3 5
	4 5 
	6 7 8
	...

there number i in the j-th row of the file indicates a nonzero of the matrix at (i, j).

It outputs a single number, namly the rank.

Compile with 

	make compile

The output is in the folder UV.

This routine is called from topology_localcopy.py.
