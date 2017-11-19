
# RA, 2017-11-17

# Plot the results of exp3b

### IMPORTS -- #

import pickle
import numpy    as np
import networkx as nx
import matplotlib.pyplot as plt

### INPUT ---- #

input_file_stats = "./OUTPUT/column-stratify-betti-3b-rand.pkl"

### OUTPUT --- #

output_file_plot = input_file_stats + "." # + eps/png

### MEAT ----- #

# LL is a list of lists
# Append zeros to each list for uniform length
def padzeros(LL) :
	maxL = max(len(L) for L in LL)
	return [(L + ([0] * (maxL - len(L)))) for L in LL]

# 
data = pickle.load(open(input_file_stats, "rb"))
print("Data keys:", list(data.keys()))

# Betti number stats (mean and stddev)
BEM = np.vstack(padzeros(data['BEM']))
BES = np.vstack(padzeros(data['BES']))

# Fraction of edges kept in the graph
FE  = np.asarray(data['FE'])

for j in range(0, BEM.shape[1]) :
	#plt.errorbar(FE, BEM[:, j], yerr=BES[:, j], fmt='.-')
	plt.plot(FE, BEM[:, j], '.-')

plt.yscale('log')
plt.xscale('log')

plt.gca().set_xlim([2e-5, 2])
#plt.xticks()
plt.gca().set_ylim([2e-1, 5e7])
plt.yticks([10**e for e in range(0, 8)])

plt.xlabel("Fraction of random egdes kept")
plt.ylabel("Average Betti numbers")

plt.legend(["b{}".format(i) for i in range(BEM.shape[1])], loc='upper left')

plt.savefig(output_file_plot + "eps")
plt.savefig(output_file_plot + "png")

plt.show()
