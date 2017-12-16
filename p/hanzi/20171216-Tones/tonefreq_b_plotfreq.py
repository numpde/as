
# RA, 2017-12-16

import numpy as np
import matplotlib.pyplot as plt

IFILE = {
	'char tone' : "OUTPUT/a/char-tones.csv"
}

OFILE = {
	'freq plot' : "OUTPUT/b/tonefreq.{ext}"
}

# Read the chararacter-tone list
# The list is sorted by frequency of the character
CT = [
	(L[0], list(int(t) for t in L[1].split(' '))) 
	for L in [
		L.rstrip().split('\t') 
		for L in open(IFILE['char tone'], "r").readlines()
	]
]

# Prepare the tone frequency table
F = np.zeros((len(CT), 5), dtype='d')
# Mark tone occurrences
for (n, (_, T)) in enumerate(CT) : F[n, T] += 1
# Count the same
F = np.cumsum(F, axis=0)

# Plot and save
plt.plot(range(len(CT)), F)
plt.xlabel("Number of most frequent characters")
plt.ylabel("Tone occurrences")
plt.legend(range(5))
plt.savefig(OFILE['freq plot'].format(ext="png"))
