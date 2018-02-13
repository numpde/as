
# RA, 2018-02-07

# Run as
#    python3 j*.py


## ================== IMPORTS :

import os
import sys
import math
import pickle
import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## ==================== INPUT :

IFILE = {
	# METABRIC expression data (log-intensity)
	'BRIC' : "OUTPUT/e_prepared/UV/bric.pkl",
	
	# TCGA expression data (FPKM)
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
}


## =================== OUTPUT :

OFILE = {
}

# Create output directories
for f in OFILE.values() : os.makedirs(os.path.dirname(f), exist_ok=True)


## ==================== PARAM :

TESTMODE = ("TEST" in sys.argv)

PARAM = {
	# Record just in case
	'testmode' : TESTMODE,
	
	# Source: https://www.biostars.org/p/77590/
	# Replaced ORC6L by ORC6 (http://www.genecards.org/cgi-bin/carddisp.pl?gene=ORC6)
	'PAM50-genes' : ("ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 NDC80 NUF2 ORC6 PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T").split(),
	
	'PAM50-types' : ["Normal", "LumA", "LumB", "Her2", "Basal"],
	
	## Extension for plots
	#'ext' : { 'pdf', 'png' },
}


## ====================== AUX :

# https://stackoverflow.com/questions/34491808/how-to-get-the-current-scripts-code-in-python
THIS = inspect.getsource(inspect.getmodule(inspect.currentframe()))


## ====================== (!) :


## ===================== WORK :


def main() :
	
	X = pickle.load(open(IFILE['TCGA'], 'rb'))['X']
	Y = pickle.load(open(IFILE['BRIC'], 'rb'))['X']
	
	Y = np.exp(np.log(2) * Y)
	
	PAM50 = PARAM['PAM50-genes']
	
	X = X.ix[ X.index.isin(PAM50) ]
	Y = Y.ix[ Y.index.isin(PAM50) ]
	
	(X, Y) = X.align(Y, join='inner', axis=0)
	
	def xscore(r) :
		s = X.ix[r.name]
		return ((r - r.mean()) / r.std()) * s.std() + s.mean()

	Y = Y.apply(xscore, axis=1)
	
	#print(X)
	#print(Y)
	
	#print(Y.min())
	#exit()
	
	
	from scipy.stats import gaussian_kde
	
	for gene in PAM50 :
		
		if (gene not in X.index) : continue
		if (gene not in Y.index) : continue
		
		x = sorted(X.ix[gene].values)
		s = np.linspace(min(x), max(x), 100)
		f = gaussian_kde(x)
		
		y = sorted(Y.ix[gene].values)
		t = np.linspace(min(y), max(y), 100)
		g = gaussian_kde(y)
		
		#plt.semilogx(s, f(s), 'r-')
		hx = plt.semilogx(x, f(x), 'r.-', markersize=3)[0]
		
		#plt.semilogx(t, g(t), 'b-')
		hy = plt.semilogx(y, g(y), 'b.-', markersize=3)[0]
		
		plt.legend([hx, hy], ["TCGA", "BRIC"])
		
		plt.show()


## ==================== ENTRY :

if (__name__ == "__main__") :
	main()
