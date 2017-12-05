
# RA, 2017-11-20

# Plot the results of exp4

### IMPORTS -- #

import pickle
import numpy    as np
import networkx as nx
import matplotlib.pyplot as plt

### INPUT ---- #

input_file_stats = "./OUTPUT/exp4a_20171201-114027.pkl"

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
#print("Data keys:", list(data.keys()))

BETTI = data['BETTI']

for (fc, betti) in BETTI.items() :
	BETTI[fc] = np.vstack(padzeros(betti))

FC = []   # Fraction of max-cliques
RUN = []  # Number of runs in the statistic
BEM = []  # Mean and ...
BES = []  #      ... stddev of Betti numbers
#
for fc in sorted(BETTI.keys()) :
	FC .append(fc)
	RUN.append(len(BETTI[fc]))
	BEM.append(np.mean(BETTI[fc], 0).tolist())
	BES.append(np.std (BETTI[fc], 0).tolist())

# Betti number stats as numpy array
BEM = np.vstack(padzeros(BEM))
BES = np.vstack(padzeros(BES))

# 
for j in range(BEM.shape[1]) :
	plt.errorbar(FC, BEM[:, j], yerr=BES[:, j], fmt='.-')
	#plt.plot(FC, BEM[:, j], '.-')

plt.yscale('log')
plt.xscale('log')

#plt.gca().set_xlim([2e-5, 2])
#plt.xticks()
#plt.gca().set_ylim([2e-1, 5e7])
#plt.yticks([10**e for e in range(0, 8)])

plt.xlabel("Fraction of max-cliques")
plt.ylabel("Average Betti numbers")

plt.legend(["b{}".format(i) for i in range(BEM.shape[1])], loc='upper left')

plt.savefig(output_file_plot + "eps")
plt.savefig(output_file_plot + "png")

plt.show()
