
Step 1
===

Read the adjacency matrix from

	../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat

construct the networkx graph and write it to 

	OUTPUT/UV/column-a-graph.pkl

Find all maximal cliques in the graph and write them to

	OUTPUT/UV/column-b-maxcliques.pkl
	
Run

	python3 step1_*.py

to launch the script.

Step 2
===

From the graph file

	OUTPUT/UV/column-a-graph.pkl

compute the graph Laplacian and write it to

	OUTPUT/UV/column-d-laplacian.mat

Run

	python3 step2_*.py

to launch the script.

Step 3
===

1
---

From the graph Laplacian in

	OUTPUT/UV/column-d-laplacian.mat

compute the first 33 eigenvectors/values and write them to

	OUTPUT/UV/column-e-lapeig-unnorm.mat

Normalize the Laplacian and write its first 33 eigenvectors/values to

	OUTPUT/UV/column-e-lapeig-normal.mat

Run **step3_1_eig** in MATLAB. See the script for details. 


References:

	o) von Luxburg, "A tutorial on spectral clustering", 2007
	o) Ng, Jordan, Weiss, "On Spectral Clustering: Analysis and an algorithm", 2001

2
---

For different number of clusters, use the Ng, Jordan, Weiss spectral clustering algorithm to cluster the graph. 
Run **step3_2_kmeans** in MATLAB.
The cluster info is in the file

	OUTPUT/step3_2_kmeans/runXX/UV/column-f-clusters.mat

The script also generates "cluster flow" images that show how data points are rearranged into different clusters when the number of clusters is increased, and histograms showing cluster sizes. The images are in the same **runXX** subfolder.

3
---


The MATLAB script **step3_3_edgestat** computes the ratio of in-cluster and out-cluster edges from the graph

	../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat

using the clustering info from

	OUTPUT/step3_2_kmeans/run*/UV/column-f-clusters.mat

for all runs, and plots the ratio in/out as a function of the number of clusters. 
The results are in the folder

	OUTPUT/step3_3_edgestat/

