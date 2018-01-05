
# RA, 2017-12-21

def heatmap(D) :
	# Based on:
	# https://medium.com/@jerilkuriakose/heat-maps-with-dendrogram-using-python-d112a34e865e
	# (accessed on 2017-12-19)

	import numpy as np
	import matplotlib.pyplot as plt
	import scipy.cluster.hierarchy as sch

	# Dendrogram that comes to the left
	fig = plt.figure(figsize=(8,8))
	# Add an axes at position rect [left, bottom, width, height]
	ax1 = fig.add_axes([0.09, 0.1, 0.1, 0.7])
	Y = sch.linkage((D), method='centroid')
	# orientation='left' is reponsible for making the 
	# dendrogram appear to the left
	Z1 = sch.dendrogram(Y, orientation='left')
	ax1.set_xticks([])
	ax1.set_yticks([])

	# top side dendogram
	ax2 = fig.add_axes([0.2, 0.81, 0.7, 0.1])
	Y = sch.linkage(np.transpose(D), method='single')
	Z2 = sch.dendrogram(Y)
	ax2.set_xticks([])
	ax2.set_yticks([])

	# main heat-map
	axmatrix = fig.add_axes([0.2, 0.1, 0.7, 0.7])
	idx1 = Z1['leaves']
	idx2 = Z2['leaves']
	D = D[idx1, :]
	D = D[:, idx2]
	# the actual heat-map
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap="YlGnBu")
	axmatrix.set_xticks([])
	axmatrix.set_yticks([])

	# xticks to the right (x-axis)
	axmatrix.set_xticks(range(D.shape[1]))
	axmatrix.set_xticklabels(idx2, minor=False)
	axmatrix.xaxis.set_label_position('bottom')
	axmatrix.xaxis.tick_bottom()

	plt.xticks(rotation=-90, fontsize=8)

	# xticks to the right (y-axis)
	axmatrix.set_yticks(range(D.shape[0]))
	axmatrix.set_yticklabels(idx1, minor=False)
	axmatrix.yaxis.set_label_position('right')
	axmatrix.yaxis.tick_right()

	# to add the color bar
	axcolor = fig.add_axes([0.01, 0.1, 0.02, 0.7])
	cbar = fig.colorbar(im, cax=axcolor)


import numpy as np
import matplotlib.pyplot as plt

x = np.random.rand(40)
D = np.zeros([20,40])
for i in range(20):
    for j in range(40):
        D[i, j] = abs(x[i] - x[j])

heatmap(D)
plt.show()
