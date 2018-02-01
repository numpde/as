
# http://zachmoshe.com/2017/04/03/pickling-keras-models.html

import types
import tempfile
import keras.models

def make_keras_picklable():
	def __getstate__(self):
		model_str = ""
		with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=True) as fd:
			keras.models.save_model(self, fd.name, overwrite=True)
			model_str = fd.read()
		d = { 'model_str': model_str }
		return d

	def __setstate__(self, state):
		with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=True) as fd:
			fd.write(state['model_str'])
			fd.flush()
			model = keras.models.load_model(fd.name)
		self.__dict__ = model.__dict__

	keras.models.Model.__getstate__ = __getstate__
	keras.models.Model.__setstate__ = __setstate__

make_keras_picklable()
