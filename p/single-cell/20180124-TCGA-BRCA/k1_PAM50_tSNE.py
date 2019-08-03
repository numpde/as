
# RA, 2019-08-02

from helpers import commons

import matplotlib.pyplot as plt
import pandas as pd
import pickle
from sklearn.manifold import TSNE
from itertools import groupby

PARAM = {
	'TCGA': "OUTPUT/e_prepared/UV/tcga.pkl",

	# https://www.biostars.org/p/77590/
	# Changed ORC6L to ORC6 (https://en.wikipedia.org/wiki/ORC6, 2019-08-02)
	'PAM50': list((
		"ACTR3B, ANLN, BAG1, BCL2, BIRC5, BLVRA, CCNB1, CCNE1, CDC20, CDC6, CDH3, CENPF, CEP55, CXXC5, EGFR, ERBB2, ESR1, EXO1, FGFR4, FOXA1, FOXC1, GPR160, GRB7, KIF2C, KRT14, KRT17, KRT5, MAPT, MDM2, MELK, MIA, MKI67, MLPH, MMP11, MYBL2, MYC, NAT1, NDC80, NUF2, ORC6, PGR, PHGDH, PTTG1, RRM2, SFRP1, SLC39A6, TMEM45B, TYMS, UBE2C, UBE2T"
	).replace(' ', '').split(',')),

	# PAM50 labels order
	'PAM50-labels' : ["Normal", "LumA", "LumB", "Her2", "Basal"],

	# Plot markers
	'PAM50-markers': ['o', '^', 'v', '*', 's'],

	# tSNE plot save-to path
	'fig': "OUTPUT/k1_PAM50_tSNE/tsne{:0>2d}.{ext}",

	# Number of tSNE repeats
	'tSNE-repeats': 6,
}


def zscore(X: pd.DataFrame, axis: int) -> pd.DataFrame:
	return X.apply((lambda s: (s - s.mean()) / (s.std() or 1)), axis=axis)


def get_data():
	TCGA_DATA = pickle.load(open(PARAM['TCGA'], 'rb'))

	# Expression table (HGNC Symbol x Sample ID)
	# Transpose to a table (Sample ID x HGNC Symbol)
	# Restrict to the PAM50 genes
	X: pd.DataFrame
	X = zscore(TCGA_DATA['X'].loc[PARAM['PAM50']].T, axis=1)

	# Check: each row sums to zero
	assert(all((abs(X.loc[i].sum()) < 1e-10) for i in X.index))

	# PAM50 assessment (Sample ID => PAM50 type)
	C: pd.Series
	C = TCGA_DATA['subtype'].set_index('aliquot_barcode', verify_integrity=True)['PAM50']

	# Restrict to the Sample ID common to both
	(X, C) = X.align(C, join='inner', axis=0)

	return (X, C)


def main():
	# Features, reference class
	(X, C) = get_data()

	# PAM50 type --> Index group
	groups = C.index.groupby(C)

	assert(set(PARAM['PAM50-labels']) == groups.keys())

	for repeat in range(PARAM['tSNE-repeats']):

		(x, y) = TSNE(n_components=2).fit_transform(X).T

		x = pd.Series(data=x, index=C.index)
		y = pd.Series(data=y, index=C.index)

		fig: plt.Figure
		ax: plt.Axes
		(fig, ax) = plt.subplots()

		colors = { pam: ("C{}".format(n)) for (n, pam) in enumerate(PARAM['PAM50-labels']) }
		markers = dict(zip(PARAM['PAM50-labels'], PARAM['PAM50-markers']))

		for pam in sorted(PARAM['PAM50-labels'], key=(lambda pam: -len(groups[pam]))):
			ax.scatter(x[groups[pam]], y[groups[pam]], c=colors[pam], label=pam, s=16, marker=markers[pam])

		lim = tuple(f(a, b) for (a, b, f) in zip(ax.get_xlim(), ax.get_ylim(), [min, max]))

		ax.set_xlim(lim)
		ax.set_ylim(lim)

		ax.legend(loc='upper left')

		ax.set_xticks([])
		ax.set_yticks([])

		ax.set_aspect(1)
		ax.grid()

		fig.savefig(
			commons.makedirs(PARAM['fig'].format(repeat, ext='png')),
			bbox_inches='tight', pad_inches=0,
			dpi=300
		)

		plt.close(fig)


if (__name__ == "__main__"):
	main()
