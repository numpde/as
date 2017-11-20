
# RA, 2017-11-07

import pickle
import numpy as np
import matplotlib.pyplot as plt

# INPUT
input_file_stats = "gnpbetti-out.pkl"

# OUTPUT
output_file_plot = "gnpbetti.eps"


# Read the statistics
data = pickle.load(open(input_file_stats, "rb"))
p = data['P']

M = np.asmatrix(data['M'])
S = np.asmatrix(data['S'])

# Number of Betti numbers
nb = max(M.nonzero()[1])

for b in range(0, nb) :
	plt.errorbar(p, M[:, b], yerr = S[:, b])


plt.legend(["b{}".format(n) for n in range(0, 10)])
plt.savefig(output_file_plot)
