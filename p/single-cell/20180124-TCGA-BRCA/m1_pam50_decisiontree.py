
# RA, 2019-08-02

# Retrain the mock PAM50 classifier on the TCGA-BRCA dataset


from helpers import commons
from commons import makedirs
from commons import logger

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

	# PAM50 label order
	'PAM50-labels' : ["Normal", "LumA", "LumB", "Her2", "Basal"],

	# Output: Regression model
	'model': makedirs("OUTPUT/m1_pam50_decisiontree/rf_model.dat"),

	# Output: Confusion matrices and accuracy score
	'confusion': makedirs("OUTPUT/m1_pam50_decisiontree/confusion{set}.txt"),
	'accuracy': makedirs("OUTPUT/m1_pam50_decisiontree/accuracy{set}.txt"),
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

	forest = RandomForestClassifier(
		n_estimators=111, max_depth=6,
		criterion='entropy', class_weight='balanced',
		random_state=1
	).fit(X0, y0)

	logger.info("Train set acc: {}".format(accuracy_score(y0, forest.predict(X0))))
	logger.info("Test set acc: {}".format(accuracy_score(y1, forest.predict(X1))))

	def save_training_report(Xi, yi, i):
		with open(PARAM['confusion'].format(set=i), 'w') as fd:
			cm = commons.confusion(yi, forest.predict(Xi), PARAM['PAM50-labels'])
			print(cm, file=fd)
		with open(PARAM['accuracy'].format(set=i), 'w') as fd:
			sc = accuracy_score(yi, forest.predict(Xi))
			print(sc, file=fd)

	save_training_report(X0, y0, 0)
	save_training_report(X1, y1, 1)

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
		open(makedirs(PARAM['model'].format()), 'wb')
	)


if (__name__ == "__main__"):
	main()
