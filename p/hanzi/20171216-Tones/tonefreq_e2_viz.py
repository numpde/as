
# https://github.com/raghakot/keras-vis/blob/master/examples/mnist/attention.ipynb

from keras.models import load_model


from scipy.misc import imsave
import numpy as np
import time
from keras.applications import vgg16
from keras import backend as K


import matplotlib.pyplot as plt


IFILE = {
	'XY' : "OUTPUT/c/UV/XY.pkl"
}

import pickle
data = pickle.load(open(IFILE['XY'], 'rb'))

(X, Y) = (data['X'], data['Y'])
#print(X.shape, Y.shape)

X = X / np.max(X)

# https://machinelearningmastery.com/handwritten-digit-recognition-using-convolutional-neural-networks-python-keras/
X = X.reshape(X.shape[0], 1, 24, 24).astype('float32')


#plt.imshow(X[2][0], cmap='gray')
#plt.show()


K.set_image_data_format('channels_first')
#print(K.image_data_format())

# dimensions of the generated pictures for each filter.
img_width = 24
img_height = 24

import sys
output_index = int(sys.argv[1])
assert(output_index in [0, 1, 2, 3, 4]), "Invalid tone: {}".format(output_index)

# util function to convert a tensor into a valid image


# 
model = load_model("OUTPUT/d/UV/model-01000.h5")
#print('Model loaded.')
#model.summary()



from vis.visualization import visualize_activation
from vis.visualization import visualize_saliency
from vis.utils import utils
from keras import activations

# Utility to search for layer index by name. 
# Alternatively we can specify this as -1 since it corresponds to the last layer.
layer_idx = utils.find_layer_idx(model, 'dense_2')
#print(model.layers[layer_idx])

### Swap softmax with linear
#model.layers[layer_idx].activation = activations.linear
#model = utils.apply_modifications(model)

#s = X[2]
#grads = visualize_saliency(model, layer_idx, filter_indices=output_index, seed_input=s, backprop_modifier='relu')
## Plot with 'jet' colormap to visualize as a heatmap.
#plt.imshow(grads, cmap='gray')
#print(grads)
#plt.show()


img = np.random.rand(*(X[2].shape))
#img = X[2]
for i in range(20) :
	img = visualize_activation(model, layer_idx, filter_indices=output_index, input_range=(0, 255), tv_weight=1, lp_norm_weight=10, seed_input=img)
#print(img.shape)
plt.imshow(img[..., 0], cmap='gray')
plt.show()



