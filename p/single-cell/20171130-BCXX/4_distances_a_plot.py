
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
	plt.clf()
	n = 0
	for (i, M) in enumerate(measurements[dims]) :
		for (j, N) in enumerate(measurements[dims]) :
			n = n + 1
			L = len(measurements[dims])
			plt.subplot(L, L, n)
			
			plt.plot(M, N, '.', markersize=1)
			
			plt.gca().set_xticks([])
			plt.gca().set_yticks([])
	
	plt.savefig(output_file_plot.format(dims=dims, m1="all", m2="all", extension="png"))

labels = []
for m in metrics :
	if (m['func'] == "ks") :
		label = "Average KS distance"
	else :
		label = m['data'] + "-" + m['func'] + " "
		if (m['info'] == "sc") : label += "(average silh. value)"
		if (m['info'] == "fr") : label += "(fraction of positive silh. values)"
	
	labels.append(label)

for dims in DIMS : 
	for (i, M) in enumerate(measurements[dims]) :
		for (j, N) in enumerate(measurements[dims]) :
			plt.clf()
			plt.plot(M, N, '.')
			plt.xlabel(labels[i])
			plt.ylabel(labels[j])
			plt.savefig(output_file_plot.format(dims=dims, m1=i, m2=j, extension="png"))
