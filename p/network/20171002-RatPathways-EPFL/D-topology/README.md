The module **topology.py** implements variants of a routine to compute the Betti numbers of graphs,
see [A. Zomorodian, Computational topology (Notes), 2009](http://www.ams.org/meetings/short-courses/zomorodian-notes.pdf). The variants differ in how the rank of the matrix of the boundary operator is computed. 

The variant **betti** uses the SVD. In the subfolder **correctness**, this it is verified against an independent implementation.

In the variant **betti_bin**, the matrix is considered as a binary matrix and the rank is computed over Z/2Z.
The result is, however, the same for simplicial complexes that arise from graphs.
This is because no two simplices share more than one face, so that the faces can be oriented consistently.
The advantage of the "binary" point of view is that only the sparsity pattern has to be stored and manipulated, not the entries of the matrix. The rank is computed via Gaussian elimination; note that, over Z/2Z, every nonzero entry is a possible pivot, but the choice of the pivot can have a dramatic impact on fill-in in the long run.

In the variant **betti_bin_cpp**, an external c++ implementation of the "binary rank" routine is invoked (see subfolder **cpp**). As of 2017-11-21, **variants/rank-9.cpp** is generally the fastest implementation, as it keeps the fill-in during the elimination relatively low (because of the lexicographic ordering of simplices produced by the python caller routine) and uses an adaptive number of computing threads for parallelism.
