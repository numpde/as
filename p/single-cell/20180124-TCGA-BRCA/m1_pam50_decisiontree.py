
# RA, 2019-08-02

# Retrain the mock PAM50 classifier on the TCGA-BRCA dataset


from helpers import commons

import pandas as pd
import json, pickle
import datetime

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score


PARAM = {
	# Input: Preprocessed TCGA-BRCA dataset with PAM50 labels
	'TCGA': "OUTPUT/e_prepared/UV/tcga.pkl",

	# List of PAM50 genes
	'PAM50': json.load(open("ORIGINALS/GeneSets/pam50/pam50.json", 'r'))['data'],

	# Output: Regression model
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

	(X, X1, y, y1) = train_test_split(X, y, test_size=0.3, random_state=1)

	forest = RandomForestClassifier(
		n_estimators=111, max_depth=6,
		criterion='entropy', class_weight='balanced',
		random_state=1
	).fit(X, y)

	M0 = commons.confusion(y, forest.predict(X), PARAM['PAM50-labels'])
	print(M0, file=open(commons.makedirs(PARAM['confusion'].format(set=0)), 'w'))
	commons.logger.info("Train set acc: {}".format(accuracy_score(y, forest.predict(X))))

	M1 = commons.confusion(y1, forest.predict(X1), PARAM['PAM50-labels'])
	print(M1, file=open(commons.makedirs(PARAM['confusion'].format(set=1)), 'w'))
	commons.logger.info("Test set acc: {}".format(accuracy_score(y1, forest.predict(X1))))

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
