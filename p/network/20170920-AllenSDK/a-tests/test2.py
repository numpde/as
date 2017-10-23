
# PYTHON 2 VERSION

from __future__ import print_function

# Path to Allen SDK
import sys
sys.path.append('../git2/')


# Need to:
# source ~/virtualenv/py2a/bin/activate
# then run
# python test2.py



from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

# The manifest file is a simple JSON file that keeps track of all of
# the data that has already been downloaded onto the hard drives.
# If you supply a relative path, it is assumed to be relative to your
# current working directory.
mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')

# open up a list of all of the experiments
all_experiments = mcc.get_experiments(dataframe=True)

# take a look at what we know about an experiment with a primary motor injection
#all_experiments.loc[122642490]





# grab the StructureTree instance
structure_tree = mcc.get_structure_tree()




from allensdk.api.queries.ontologies_api import OntologiesApi

oapi = OntologiesApi()

# get the ids of all the structure sets in the tree
structure_set_ids = structure_tree.get_structure_sets()

# query the API for information on those structure sets
#print(oapi.get_structure_sets(structure_set_ids))

# List of Primary injection structures for BDA/A... 	114512892 	Mouse Connectivity - BDA/AAV Primary Injection...
# List of primary AND secondary injection struct... 	112905813 	Mouse Connectivity - BDA/AAV All Injection Str...
# List of structures for ABA Fine Structure Search 	10 	ABA - Fine Structure Search
# List of structures used for the Connectivity p... 	167587189 	Mouse Connectivity - Summary
# List of primary AND secondary injection struct... 	112905828 	Mouse Connectivity - Projection All Injection ...
# All mouse visual areas with layers 	396673091 	Mouse Cell Types - Structures
# None 	514166994 	CAM targeted structure set
# List of structures for ABA Differential Search 	12 	ABA - Differential Search
# List of structures representing a coarse level... 	2 	Mouse - Coarse
# List of valid structures for projection target... 	184527634 	Mouse Connectivity - Target Search
# List of structures representing a areal level ... 	3 	Mouse - Areas
# List of Primary injection structures for Proje... 	114512891 	Mouse Connectivity - Projection Primary Inject...





structures = structure_tree.get_structures_by_set_id([112905813])
#print(structures)
for S in structures : 
	if ('cortex' in S['name']) :
		 print(S['id'], S['name'])

# 184 Frontal pole, cerebral cortex
# 315 Isocortex
# 528 Cerebellar cortex
# 688 Cerebral cortex
# 856 Thalamus, polymodal association cortex related
# 864 Thalamus, sensory-motor cortex related






# fetch the experiments that have injections in the isocortex of cre-positive mice
isocortex = structure_tree.get_structures_by_name(['Isocortex'])[0]
print("Isocortex ID:", isocortex['id'])

cre_cortical_experiments = mcc.get_experiments(cre=True, injection_structure_ids=[184])
print(cre_cortical_experiments[0])






# PROJECTION MATRIX


#import matplotlib.pyplot as plt

for S in structures : 
	isid = S['id'] # Injection structure ID
	#isid = structure_tree.get_structures_by_acronym(['VISp'])[0]['id']
	
	print("Injection structure ID:", isid)

	experiments = mcc.get_experiments(cre = False, injection_structure_ids = [isid])
	experiment_ids = [e['id'] for e in experiments]

	print("Experiment IDs:", experiment_ids)

	ctx_children = structure_tree.child_ids( [isocortex['id']] )[0]
	print("ctx_children =", ctx_children)
	pm = mcc.get_projection_matrix(experiment_ids = experiment_ids, projection_structure_ids = ctx_children,
                               hemisphere_ids = [2], # right hemisphere, ipsilateral
                               parameter = 'projection_density')
	
	print(pm)

	'''
	row_labels = pm['rows'] # these are just experiment ids
	column_labels = [ c['label'] for c in pm['columns'] ] 
	matrix = pm['matrix']

	fig, ax = plt.subplots(figsize=(15,15))
	heatmap = ax.pcolor(matrix, cmap=plt.cm.afmhot)

	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(matrix.shape[1])+0.5, minor=False)
	ax.set_yticks(np.arange(matrix.shape[0])+0.5, minor=False)

	ax.set_xlim([0, matrix.shape[1]])
	ax.set_ylim([0, matrix.shape[0]])          

	# want a more natural, table-like display
	ax.invert_yaxis()
	ax.xaxis.tick_top()

	ax.set_xticklabels(column_labels, minor=False)
	ax.set_yticklabels(row_labels, minor=False)
	plt.figure()
	plt.ion()
	plt.show()
	'''
	

exit(0)





# find wild-type injections into primary visual area
#visp = structure_tree.structures_get_by_acronym(['VISp'])[0]
visp_experiments = mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[visp['id']])





import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
#matplotlib inline

visp_experiment_ids = [ e['id'] for e in visp_experiments ]
ctx_children = structure_tree.child_ids( [isocortex['id']] )[0]

pm = mcc.get_projection_matrix(experiment_ids = visp_experiment_ids, 
                               projection_structure_ids = ctx_children,
                               hemisphere_ids= [2], # right hemisphere, ipsilateral
                               parameter = 'projection_density')

row_labels = pm['rows'] # these are just experiment ids
column_labels = [ c['label'] for c in pm['columns'] ] 
matrix = pm['matrix']

fig, ax = plt.subplots(figsize=(15,15))
heatmap = ax.pcolor(matrix, cmap=plt.cm.afmhot)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(matrix.shape[1])+0.5, minor=False)
ax.set_yticks(np.arange(matrix.shape[0])+0.5, minor=False)

ax.set_xlim([0, matrix.shape[1]])
ax.set_ylim([0, matrix.shape[0]])          

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(column_labels, minor=False)
ax.set_yticklabels(row_labels, minor=False)
plt.show()


