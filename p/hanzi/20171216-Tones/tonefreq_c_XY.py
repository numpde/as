
# RA, 2017-12-16

IFILE = {
	'char tone' : "OUTPUT/a/char-tones.csv",
	'font'      : "ORIGINALS/UV/NotoSansCJK-Regular.ttc"
}

OFILE = {
	'XY' : "OUTPUT/c/UV/XY.pkl"
}

PARAM = {
	'image size' : 24
}

# Read the chararacter-tone list
# The list is sorted by frequency of the character

CT = [
	(L[0], list(int(t) for t in L[1].split(' '))) 
	for L in [
		L.rstrip().split('\t') 
		for L in open(IFILE['char tone'], "r").readlines()
	]
]

# Generate images of characters 

from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt
import numpy as np

w = PARAM['image size']
X = np.zeros((len(CT), w, w))
Y = np.zeros((len(CT), 5))

for (n, (s, T)) in enumerate(CT) :
	
	filename = '{0:04d}'.format(n)
	
	font = ImageFont.truetype(IFILE['font'], w)
	(_, h) = font.getsize(s)
	
	(_, h) = font.getsize(s)
	img = Image.new('L', (w, w), 0xffffff)
	ImageDraw.Draw(img).text( (0, w-h), s, font = font, fill = 0x000000 )
	
	I = np.array(img)
	
	#plt.ion()
	#plt.imshow(I, cmap='grey')
	#plt.show()
	
	X[n, :, :] = I
	
	T = list(set(T))
	Y[n, T] = 1 / len(T)


# Save images of characters and their tones

import pickle
pickle.dump({ 'X' : X, 'Y' : Y }, open(OFILE['XY'], 'wb'))

