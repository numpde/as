The script **extract.py** reads 
the neural connectome from
the Blue Brain Project file **pathways_mc0_Column.h5** (located under **ORIGINALS**),
and
writes the graph as an adjacency matrix and 
the 3d locations of neurons in 3d into 
the MATLAB file **pathways_mc0_Column.h5.mat** in the **OUTPUT** folder.

In addition it writes the incidence list to **ratcolumn_data.txt** with some metainfo in **ratcolumn_head.txt**.

Run

	make x

to launch the script.
