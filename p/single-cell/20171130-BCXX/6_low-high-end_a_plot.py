
# RA, 2017-12-12

# Empirical frequency of GO terms by clustering index
# over random samples of genes.

## ================== IMPORTS :

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
from scipy.stats     import gaussian_kde

## ==================== INPUT :

IFILE = {
	"Results" : "OUTPUT/6_low-high-end_a/UV/results.pkl",
	"GO info" : "ORIGINALS/go/go-summary.csv",
}

## =================== OUTPUT :

OFILE = {
	"Frequency" : "OUTPUT/6_low-high-end_a/go-freq_dims={dims}_filt={filt}.{ext}"
}

## =================== PARAMS :

pass

## ==================== PREPA :

# Check existence of input files
for f in IFILE.values() :
	assert(os.path.isfile(f)), "File {} not found.".format(f)

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)
	
## ===================== WORK :


def main() :
	
	# [ LOAD GO ID ANNOTATION ]
	
	GO2T = dict(
		tuple(L.rstrip().split('\t')[:2])
		for L in open(IFILE["GO info"], 'r').readlines()[1:]
	)

	# [ LOAD RESULTS ]

	GO2I = pickle.load(open(IFILE["Results"], "rb"))['GO2I']
	
	# [ PLOT ]
	
	for dims in sorted(GO2I.keys()) :
		
		# Remove possible NaN and Inf
		GO2I[dims] = { go : [i for i in I if np.isfinite(i)] for (go, I) in GO2I[dims].items() }

		# 
		lenI = list(reversed(sorted(len(I) for I in GO2I[dims].values())))
		
		for min_len_rank in [4, 16, 64] :
			
			print("dims = {}, rank = {}".format(dims, min_len_rank))
			
			GO_I = []
			
			for (go, I) in GO2I[dims].items() :
				
				if (len(I) == 0) : continue
				if (min(I) == max(I)) : continue
				if (len(I) < lenI[min_len_rank]) : continue
				
				GO_I.append( (go, I) )
			
			# Sort by clustering index, LOWEST MEDIAN FIRST
			GO_I = sorted(((go, I) for (go, I) in GO_I), key=(lambda go_I : np.median(go_I[1])))
			
			plt.clf()
			
			for (go, I) in GO_I :
				
				t = np.linspace(min(I), max(I), 100)
				f = gaussian_kde(I)(t)
				
				plt.plot(t, f)
			
			plt.legend([GO2T.get(go, go) for (go, I) in GO_I], loc='upper left', prop={'size': 6})
			plt.xlabel("Clustering index")
			plt.ylabel("Relative empirical frequency")
			plt.savefig(OFILE["Frequency"].format(dims=dims, filt=min_len_rank, ext="png"))

if (__name__ == "__main__") :
	main()
