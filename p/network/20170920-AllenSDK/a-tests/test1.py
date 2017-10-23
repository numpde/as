import sys
sys.path.append('../git/')

# https://github.com/AllenInstitute/AllenSDK/blob/master/doc_template/examples/nb/mouse_connectivity.html

global xrange
def xrange(r) : return range(r)
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

# The manifest file is a simple JSON file that keeps track of all of
# the data that has already been downloaded onto the hard drives.
# If you supply a relative path, it is assumed to be relative to your
# current working directory.
mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')

# open up a list of all of the experiments
all_experiments = mcc.get_experiments(dataframe=True)
print("{} total experiments".format(len(all_experiments)))

# take a look at what we know about an experiment with a primary motor injection
print(all_experiments.loc[122642490])

# grab the StructureTree instance
structure_tree = mcc.get_structure_tree()

print(structure_tree)


from allensdk.api.queries.ontologies_api import OntologiesApi

oapi = OntologiesApi()

# get the ids of all the structure sets in the tree
structure_set_ids = structure_tree.get_structure_sets()

# query the API for information on those structure sets
#print(oapi.get_structure_sets(structure_set_ids))




#
#structures = structure_tree.get_structures_by_set_id([167587189])
#print(structures)








# Projection grid data volume

experiment_id = 181599674

# projection density: number of projecting pixels / voxel volume
pd, pd_info = mcc.get_projection_density(experiment_id)

# injection density: number of projecting pixels in injection site / voxel volume
ind, ind_info = mcc.get_injection_density(experiment_id)

# injection fraction: number of pixels in injection site / voxel volume
inf, inf_info = mcc.get_injection_fraction(experiment_id)

# data mask:
# binary mask indicating which voxels contain valid data
dm, dm_info = mcc.get_data_mask(experiment_id)

template, template_info = mcc.get_template_volume()
annot, annot_info = mcc.get_annotation_volume()

print(pd_info)
print(pd.shape, template.shape, annot.shape)

