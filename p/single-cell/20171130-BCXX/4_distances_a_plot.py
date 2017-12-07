
# RA, 2017-12-07

## ================== IMPORTS :

import pickle
import matplotlib.pyplot as plt

## ==================== INPUT :

input_file_measurements = "OUTPUT/4_distances_a/measurements.pkl"

## =================== OUTPUT :

output_file_plot = "OUTPUT/4_distances_a/dims={dims}_m1={m1}_m2={m2}.{extension}"

## =================== PARAMS :

pass

## ===================== WORK :

data = pickle.load(open(input_file_measurements, 'rb'))

DIMS = data['DIMS']
metrics = data['metrics']
measurements = data['measurements']

for dims in DIMS : 
	for (i, M) in enumerate(measurements[dims]) :
		for (j, N) in enumerate(measurements[dims]) :
			plt.clf()
			plt.plot(M, N, '.')
			plt.xlabel(" / ".join(sorted(metrics[i].values())))
			plt.ylabel(" / ".join(sorted(metrics[j].values())))
			plt.savefig(output_file_plot.format(dims=dims, m1=i, m2=j, extension="png"))
