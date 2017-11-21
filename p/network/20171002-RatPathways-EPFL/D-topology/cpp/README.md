The c++ routine **variants/rank-X.cpp** computes the Z/2Z rank of a binary matrix, cf. [this](https://gist.github.com/numpde/9584779ad235c6ee19be7a6bb87e8af5) gist.

It reads a filename from stdin.

The file is expected to be in the format:

	1 3 5
	4 5 
	6 7 8
	...

there number **i** in the **j**-th row of the file indicates a nonzero of the matrix at **(i, j)**.

It outputs a single number, namely the rank.

Compile with 

	make compile

The executable **rank** is in the folder UV.

Test with

	make test


This routine is called from **topology.py** or **topology_localcopy.py**.
