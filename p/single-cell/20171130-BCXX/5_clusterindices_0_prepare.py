
# RA, 2017-12-12

## ================== IMPORTS :

import re
import math
import pickle
import random
import inspect
import numpy as np

from scipy           import stats
from itertools       import chain
from multiprocessing import cpu_count
from progressbar     import ProgressBar as Progress
from joblib          import Parallel, delayed

class BC :
	
	def __init__(self) :

		## ==================== INPUT :

		input_file_BC = "OUTPUT/0_select/UV/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt-selected.pkl"

		## =================== OUTPUT :

		pass

		## =================== PARAMS :

		# Number of parallel computing processes
		num_procs = min(math.ceil(cpu_count() / 2), 12)

		# Number of dots per graph
		dots = 1000

		# Number of dimensions to select for each sample (each dot)
		DIMS = [10, 100, 1000, 10000, 20000]

		# 
		random_seed = 0
		
		# Strange samples to remove
		remove_samples = { "BC07LN_20" }

		## ===================== WORK :

		# Load the BC data
		self.data = pickle.load(open(input_file_BC, "rb"))
		#print(data.keys())

		# Expression matrix
		self.X = self.data['X']

		# Labels for axis/dimension of data
		(self.axis_smpl, self.axis_gene) = (self.data['axis_smpl'], self.data['axis_gene'])

		# Labels of samples of the form BCXX[LN][_Re]_XX 
		self.sample_labels = self.data['header']

		# Remove strange samples
		for h in remove_samples :
			self.X = np.delete(self.X, self.sample_labels.index(h), axis=self.axis_smpl)
			self.sample_labels.remove(h)

		# Number of samples / genes in the expression matrix
		(self.n_samples, self.n_genes) = (self.X.shape[self.axis_smpl], self.X.shape[self.axis_gene])

		# Collect numbers of samples of the form "g"_XX by group
		self.G2S = self.data['B2SH'] # By batch
		
		self.groups = list(sorted(self.G2S.keys()))

		# Split the expression matrix by groups
		self.G2X = {
			g : np.take(self.X, self.G2S[g], axis=self.axis_smpl) 
			for g in self.G2S.keys()
		}

		self.S = list(G2S.values())

		self.Z = stats.mstats.zscore(self.X, axis=self.axis_smpl)


if (__name__ == "__main__") :
	#main()
	pass


