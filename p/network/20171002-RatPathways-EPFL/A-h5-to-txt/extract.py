import h5py
import numpy as np

# Blue Brain Project file to read
filename = "pathways_mc0_Column.h5"

f = h5py.File(filename, "r")
fc = f["connectivty"]
fp = f["populations"]

#for n in fc : print(n)
#for n in fp : print(n)

print("# Connections from {}".format(filename))

# 1. Construct the graph adjacency matrix

M = []
for (k0, v0) in fc.items() :
	mm = []
	for (k1, v1) in v0.items() :
		m = v1["cMat"]
		mm.append(m)
		#print("# Connections from {} to {}. Shape: {}".format(k0, k1, m.shape))
		del m

	M.append(np.hstack(mm))
	del mm

# This is the adjacency matrix (i -> j)
M = np.vstack(M)

print("# Adjacency matrix size: {}".format(M.shape))
#plt.ion(); plt.imshow(M); plt.show(); input("# Press enter...")


# 2a. Save in Matlab format

from scipy.sparse import csc_matrix
from scipy.io import savemat
savemat(filename + ".mat", {'M' : csc_matrix(M)}, do_compression=True)

# 2b. Output to stdout in txt format

(I, J) = np.nonzero(M)
IJ = list(zip(list(I), list(J)))

print("# Number of connections: {}".format(len(IJ)))

base = 0

print("# Array base / first node ID: {}".format(base))
print("# FromNodeId\tToNodeId")

for (i, j) in IJ :
	print("{}\t{}".format(base+i, base+j))


