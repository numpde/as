
# RA, 2019-08-02

# Retrain the mock PAM50 classifier on the TCGA-BRCA dataset


from helpers import commons

import pandas as pd
import pickle
import datetime

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics


PARAM = {
	# Preprocessed TCGA-BRCA dataset with PAM50 labels
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

	# Confusion matrices
	'confusion': "OUTPUT/m1_pam50_decisiontree/confusion{set}.txt",

	# Save regression model to disk
	'model': "OUTPUT/m1_pam50_decisiontree/rf_model.dat",
}



def get_data():
	TCGA_DATA = pickle.load(open(PARAM['TCGA'], 'rb'))

	# Expression table (HGNC Symbol x Sample ID)
	# Transpose to a table (Sample ID x HGNC Symbol)
	# Normalize to z-score sample-wise
	# Restrict features to only the PAM50 genes
	X: pd.DataFrame
	X = commons.zscore(TCGA_DATA['X'].T, axis=1)[PARAM['PAM50']]

	# # Check: each row sums to zero
	# assert(all((abs(row.sum()) < 1e-10) for (i, row) in X.iterrows()))

	# PAM50 assessment (Sample ID => PAM50 type)
	C: pd.Series
	C = TCGA_DATA['subtype'].set_index('aliquot_barcode', verify_integrity=True)['PAM50']

	# Restrict to the Sample ID common to both
	(X, C) = X.align(C, join='inner', axis=0)

	return (X, C)


def train(X, y):

	(X0, X1, y0, y1) = train_test_split(X, y, test_size=0.3, random_state=1)

	forest: RandomForestClassifier
	forest = RandomForestClassifier(
		n_estimators=111, max_depth=6,
		criterion='entropy', class_weight='balanced',
		random_state=1
	).fit(X0, y0)

	M0 = commons.confusion(y0, forest.predict(X0), PARAM['PAM50-labels'])
	print(M0, file=open(commons.makedirs(PARAM['confusion'].format(set=0)), 'w'))
	commons.logger.info("Train set accuracy: {}".format(metrics.accuracy_score(y0, forest.predict(X0))))

	M1 = commons.confusion(y1, forest.predict(X1), PARAM['PAM50-labels'])
	print(M1, file=open(commons.makedirs(PARAM['confusion'].format(set=1)), 'w'))
	commons.logger.info("Test set accuracy: {}".format(metrics.accuracy_score(y1, forest.predict(X1))))

	return forest


def main():
	# Features, reference class
	(X, C) = get_data()

	# Train the classifier
	model = train(X, C)

	pickle.dump(
		{
			'model': model,
			'features': X.columns,
			'data': {'X': X, 'C': C},

			'PARAM': PARAM,
			'script': commons.this_module_body(),
			'timestamp': datetime.datetime.utcnow(),
		},
		open(commons.makedirs(PARAM['model'].format()), 'wb')
	)


if (__name__ == "__main__"):
	main()
