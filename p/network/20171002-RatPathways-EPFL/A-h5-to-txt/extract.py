# RA, 20171023

# Input:
# Blue Brain Project file to read, h5 format
filename_input = "./ORIGINALS/UV/pathways_mc0_Column.h5"

# Output:
output_dir = "./OUTPUT"
output_filename_mat      = output_dir + "/UV/pathways_mc0_Column.h5.mat"
output_filename_rat_head = output_dir + "/UV/ratcolumn_head.txt"
output_filename_rat_data = output_dir + "/UV/ratcolumn_data.txt"

import h5py
import numpy as np

f = h5py.File(filename_input, "r")
fc = f["connectivty"]
fp = f["populations"]

#for n in fc : print(n)
#for n in fp : print(n)

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

#plt.ion(); plt.imshow(M); plt.show(); input("# Press enter...")


# 2a. Save in Matlab format

from scipy.sparse import csc_matrix
from scipy.io import savemat
savemat(output_filename_mat, {'M' : csc_matrix(M)}, do_compression=True)

# 2b. Output in txt format

(I, J) = np.nonzero(M)
IJ = list(zip(list(I), list(J)))

base = 0

with open(output_filename_rat_head, "w") as f :
    print("# Connections from {}".format(filename_input), file=f)
    print("# Adjacency matrix size: {}".format(M.shape), file=f)
    print("# Number of connections: {}".format(len(IJ)), file=f)
    print("# Array base / first node ID: {}".format(base), file=f)
    print("# FromNodeId\tToNodeId", file=f)

with open(output_filename_rat_data, "w") as f :
    for (i, j) in IJ :
            print("{}\t{}".format(base+i, base+j), file=f)


