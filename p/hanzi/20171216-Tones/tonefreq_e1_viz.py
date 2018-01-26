
# https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html
# https://github.com/keras-team/keras/issues/431

from keras.models import load_model


from scipy.misc import imsave
import numpy as np
import time
from keras.applications import vgg16
from keras import backend as K

# dimensions of the generated pictures for each filter.
img_width = 24
img_height = 24

# the name of the layer we want to visualize
# (see model definition at keras/applications/vgg16.py)
layer_name = 'conv2d_4'

import sys
output_index = int(sys.argv[1])
assert(output_index in [0, 1, 2, 3, 4]), "Invalid tone: {}".format(output_index)

# util function to convert a tensor into a valid image


def deprocess_image(x):
	# normalize tensor: center on 0., ensure std is 0.1
	x -= x.mean()
	x /= (x.std() + K.epsilon())
	x *= 0.1

	# clip to [0, 1]
	x += 0.5
	x = np.clip(x, 0, 1)

	# convert to RGB array
	x *= 255
	if K.image_data_format() == 'channels_first':
		x = x.transpose((1, 2, 0))
	x = np.clip(x, 0, 255).astype('uint8')
	return x

# 
model = load_model("OUTPUT/d/UV/model-01000.h5")

print('Model loaded.')

model.summary()

# this is the placeholder for the input images
input_img = model.input

# get the symbolic outputs of each "key" layer (we gave them unique names).
layer_dict = dict([(layer.name, layer) for layer in model.layers[1:]])


def normalize(x):
	# utility function to normalize a tensor by its L2 norm
	return x / (K.sqrt(K.mean(K.square(x))) + K.epsilon())


loss = K.mean(model.output[:, output_index])

## we build a loss function that maximizes the activation
## of the nth filter of the layer considered
#layer_output = layer_dict[layer_name].output

#if (filter_index >= layer_output.shape[1]) : break

#if K.image_data_format() == 'channels_first' or True:
	#loss = K.mean(layer_output[:, filter_index, :, :])
#else:
	#loss = K.mean(layer_output[:, :, :, filter_index])

# we compute the gradient of the input picture wrt this loss
grads = K.gradients(loss, input_img)[0]

# normalization trick: we normalize the gradient
grads = normalize(grads)

# this function returns the loss and grads given the input picture
iterate = K.function([K.learning_phase()] + [input_img], [loss, grads])

kept_filters = []
for filter_index in range(1000):
	print('Processing filter %d' % filter_index)
	start_time = time.time()

	

	# step size for gradient ascent
	step = 0.01

	# we start from a gray image with some random noise
	if K.image_data_format() == 'channels_first' or True:
		input_img_data = np.random.random((1, 1, img_width, img_height))
	else:
		input_img_data = np.random.random((1, img_width, img_height, 1))
	input_img_data = ((input_img_data - 0.5) * 20 + 128) / 100

	# we run gradient ascent for a few steps
	for i in range(10):
		loss_value, grads_value = iterate([0] + [input_img_data])
		input_img_data += grads_value * step
		
		if loss_value <= 0.:
			# some filters get stuck to 0, we can skip them
			break
	print('Loss value:', loss_value)

	# decode the resulting input image
	if loss_value > 0:
		img = deprocess_image(input_img_data[0])
		kept_filters.append((img, loss_value))
	end_time = time.time()
	print('Filter %d processed in %ds' % (filter_index, end_time - start_time))

# we will stich the best n*n filters on a n x n grid.
n = 10

# the filters that have the highest loss are assumed to be better-looking.
# we will only keep the top n*n filters.
kept_filters.sort(key=lambda x: x[1], reverse=True)
kept_filters = kept_filters[:n * n]

# build a black picture with enough space 
margin = 2
width = n * img_width + (n - 1) * margin
height = n * img_height + (n - 1) * margin
stitched_filters = np.zeros((width, height))

# fill the picture with our saved filters
for i in range(n):
	for j in range(n):
		if (len(kept_filters) <= (i*n + j)) : continue
		img, loss = kept_filters[i * n + j]
		stitched_filters[(img_width + margin) * i: (img_width + margin) * i + img_width,
						(img_height + margin) * j: (img_height + margin) * j + img_height] = img.reshape((img_width, img_height))


imsave('OUTPUT/e1/excites_tone_{}.png'.format(output_index), stitched_filters)
