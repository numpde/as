
# RA, 2017-12-16

IFILE = {
	'char freq' : "ORIGINALS/CharFreq-Modern.csv"
}

OFILE = {
	'char tone' : "OUTPUT/a/char-tones.csv"
}

# Read lines, split by TAB
LL = [L.rstrip().split('\t') for L in open(IFILE['char freq'], "r").readlines()[1:]]

# Character-Pronunciation list
CP = list((L[1], L[4]) for L in LL if (len(L) > 4))

# Convert something like "de/di2/di4" to a list of tones [0, 2, 4]
def p2t(T) :
	return [int('0' + ''.join(filter(str.isdigit, t))) for t in T.split('/')]

# Convert each pronunciation entry to a list of tones
CT = [(c, p2t(p)) for (c, p) in CP]

# Write the character-tones list to file
with open(OFILE['char tone'], "w") as f :
	for (c, T) in CT :
		print(c + '\t' + ' '.join(str(t) for t in T), file=f)
