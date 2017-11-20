
# RA, 2017-11-13

# Plot the results of exp2a

### IMPORTS -- #

import pickle
import numpy    as np
import networkx as nx
import matplotlib.pyplot as plt

### INPUT ---- #

input_file_stats = "./OUTPUT/column-stratify-stats-2a-dist.pkl"

### OUTPUT --- #

output_file_plot = input_file_stats + "." # + eps/png

### MEAT ----- #

# LL is a list of lists
# Append zeros to each list for uniform length
def padzeros(LL) :
	maxL = max(len(L) for L in LL)
	return [(L + ([0] * (maxL - len(L)))) for L in LL]

# data has keys ['maxd', 'P', 'FE', 'NC']
data = pickle.load(open(input_file_stats, "rb"))

# We omit 0-cliques
NC = np.vstack(padzeros(data['NC']))[:, 1:]


FE  = np.asarray(data['FE'])


plt.plot(FE, NC, '.-')


plt.yscale('log')
plt.xscale('log')

plt.gca().set_xlim([2e-5, 2])
#plt.xticks()
plt.gca().set_ylim([2e-1, 2e8])
plt.yticks([10**e for e in range(0, 9)])

plt.xlabel("Fraction of shortest egdes kept")
plt.ylabel("Number of cliques")

plt.legend(["{}-cliques".format(1+i) for i in range(NC.shape[1])], loc='upper left')

plt.savefig(output_file_plot + "eps")
plt.savefig(output_file_plot + "png")

plt.show()
