
# RA, 2017-11-07

import pickle
import numpy as np
import matplotlib.pyplot as plt

# INPUT
input_file_stats = "./OUTPUT/gnpbetti-out_n={}.pkl"

# OUTPUT
output_file_plot = "./OUTPUT/gnpbetti_n={}"

for n in [10, 20, 30, 40, 60, 70, 80] :
	
	# Read the statistics
	data = pickle.load(open(input_file_stats.format(n), "rb"))
	p = np.asmatrix(data['P']).reshape(-1, 1)

	M = np.asmatrix(data['M'])
	S = np.asmatrix(data['S'])

	# Number of Betti numbers
	nb = max(M.nonzero()[1])

	plt.clf()
	
	for b in range(0, nb) :
		plt.errorbar(p, M[:, b], yerr = [M[:, b] - np.maximum(M[:, b] - S[:, b], 0), S[:, b]])

	plt.xlabel("Edge-inclusion probability p for a graph on {} nodes".format(n))
	plt.ylabel("Betti number")

	plt.legend(["0th", "1st", "2nd", "3rd"] + ["{}th".format(k) for k in range(4, 10)])
	plt.savefig(output_file_plot.format(n) + ".eps")
	plt.savefig(output_file_plot.format(n) + ".png")
