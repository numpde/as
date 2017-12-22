
# RA, 2017-12-20

## ================== IMPORTS :

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
from scipy.stats     import gaussian_kde

## ==================== INPUT :

IFILE = {
	"Results" : "OUTPUT/6_low-high-end_b/UV/results.pkl",
	"GO info" : "ORIGINALS/go/go-summary.csv",
	"GO -> ENSG" : "OUTPUT/0_e2go/go2e.txt",
}

## =================== OUTPUT :

OFILE = {
	"Frequency" : "OUTPUT/6_low-high-end_b/go-freq_dims={dims}.{ext}"
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
	
	# GO2T : GO ID --> Text label
	GO2T = dict(
		tuple(L.rstrip().split('\t')[:2])
		for L in open(IFILE["GO info"], 'r').readlines()[1:]
	)
	
	# GO2N : GO ID --> Number of genes
	GO2N = {

		go_E[0] : len(go_E[1:])

		for go_E in [
			L.rstrip().split('\t') 
			for L in open(IFILE["GO -> ENSG"], 'r')
		]
	}

	# [ LOAD RESULTS ]

	GO2I = pickle.load(open(IFILE["Results"], "rb"))['GO2I']
	
	# [ PLOT ]
	
	for dims in sorted(GO2I.keys()) :
		
		# Extract x & y for plotting
		
		GO_I = []
		
		for (go, I) in GO2I[dims].items() :
			# Remove possible NaN and Inf
			I = [i for i in I if np.isfinite(i)]
			
			if (len(I) == 0) : continue
			if (min(I) == max(I)) : continue
			
			GO_I.append( (go, I) )
		
		# Sort by clustering index, LOWEST MEDIAN FIRST
		GO_I = sorted(((go, I) for (go, I) in GO_I), key=(lambda go_I : np.median(go_I[1])))
		
		# Take the most interesting ones
		GO_I = GO_I[0:200]
		
		# Plot
		
		plt.clf()
		
		colors = plt.get_cmap('hsv')(np.linspace(0.1, 1.0, len(GO_I))).tolist()
		
		for (n, (go, I)) in enumerate(GO_I) :
			
			t = np.linspace(min(I), max(I), 100)
			f = gaussian_kde(I)
			
			plt.plot(t, f(t), '-', color=colors.pop(), zorder=(-5*n))
		
		
		# Legend labels
		L = [(GO2T.get(go, go) + " ({})".format(GO2N[go])) for (go, I) in GO_I]
		# and font size
		legend_font_size = 4
		
		plt.legend(L[0:40], loc='upper left', prop={'size': legend_font_size})
		plt.xlabel("Clustering index")
		plt.ylabel("Relative empirical frequency")
		
		plt.savefig(OFILE["Frequency"].format(dims=dims, ext="png"))
		plt.savefig(OFILE["Frequency"].format(dims=dims, ext="eps"))

if (__name__ == "__main__") :
	main()
