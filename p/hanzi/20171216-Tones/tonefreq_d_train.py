
# RA, 2017-12-16

IFILE = {
	'XY' : "OUTPUT/c/UV/XY.pkl"
}

OFILE = {
	'model' : "OUTPUT/d/UV/model-{epoch:05d}.h5"
}

import numpy as np

import pickle
data = pickle.load(open(IFILE['XY'], 'rb'))

(X, Y) = (data['X'], data['Y'])
#print(X.shape, Y.shape)

X = X / np.max(X)

# https://machinelearningmastery.com/handwritten-digit-recognition-using-convolutional-neural-networks-python-keras/
X = X.reshape(X.shape[0], 1, 24, 24).astype('float32')

num_classes = Y.shape[1]

from keras.models import Sequential
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten, Activation
from keras        import regularizers

import keras.optimizers

# https://stackoverflow.com/questions/39547279/loading-weights-in-th-format-when-keras-is-set-to-tf-format
from keras import backend as keras_backend
keras_backend.set_image_dim_ordering('th')

# The model architecture is from
# https://github.com/keras-team/keras/blob/master/examples/cifar10_cnn.py
model = Sequential()
model.add(Conv2D(32, (3, 3), padding='same', input_shape=X.shape[1:]))
model.add(Activation('relu'))
model.add(Conv2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Conv2D(64, (3, 3), padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(64, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(512))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes))
model.add(Activation('softmax'))

opt = keras.optimizers.rmsprop(lr=0.0001, decay=1e-6)

model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

# https://keras.io/callbacks/#example-model-checkpoints
mc = keras.callbacks.ModelCheckpoint(OFILE['model'], monitor='loss', verbose=0, save_best_only=False, save_weights_only=False, mode='auto', period=1)

model.fit(X, Y, epochs=1000, callbacks=[mc])
