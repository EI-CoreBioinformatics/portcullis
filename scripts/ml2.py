#!/usr/bin/env python3

import os
import argparse
import bed12
from sklearn.cross_validation import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import classification_report

import numpy as np
import matplotlib.pyplot as plt



def main():
	parser = argparse.ArgumentParser("Script to build a random forest decision tree")
	parser.add_argument("pos", help="The tab file produce by portcullis")
	parser.add_argument("neg", help="The tab file produce by portcullis")
	parser.add_argument("input", help="The tab file produce by portcullis")
	#parser.add_argument("input2", help="The tab file produce by portcullis")
	parser.add_argument("-t", "--threads", type=int, default="1", help="The number of threads to use")
	parser.add_argument("--test", help="Test the classifier against this file")
	parser.add_argument("-o", "--output", required=True, help="The output prefix")
	args = parser.parse_args()

	# X should contain a matrix of features derived from the portcullis tab file
	# y should contain the labels (0 not a valid junction, 1 a valid junction).  Confirmed with the reference.

	pos = bed12.loadbed(args.pos, False, False)
	neg = bed12.loadbed(args.neg, False, False)


	data_X = []
	train_Y = []
	data = []
	test_Y = []
	#ref = open(args.input2)
	with open(args.input) as f:
		# Skip header
		f.readline()

		for line in f:
			parts = line.strip().split(sep="\t")
			if len(parts) > 1:
				test_Y.append(int(parts[8]))
				#test_Y.append(int(ref.readline()))
				raw = int(parts[14])
				#rel = float(parts[9])
				#rel2raw = float(parts[10])
				maxmmes = float(parts[11])

				#data.append(parts[9:-1])
				data.append(maxmmes)
				b = bed12.BedEntry(False)
				b.chrom = parts[2]
				b.thick_start = int(parts[4])
				b.thick_end = int(parts[5]) + 1
				b.strand = parts[12]
				if b in pos:
					#data_X.append(parts[9:-1])
					data_X.append(maxmmes)
					train_Y.append(1)
				if b in neg:
					#data_X.append(parts[9:-1])
					data_X.append(maxmmes)
					train_Y.append(0)


	train_X = np.array(data_X, dtype='|S4').astype(np.float).reshape(-1, 1)
	test_X = np.array(data, dtype='|S4').astype(np.float).reshape(-1, 1)
	#train_X = np.array(data_X, dtype='|S4').astype(np.float)
	#test_X = np.array(data, dtype='|S4').astype(np.float)


	print("Logistic regression with L1 regularisation")
	#L1 regularized logistic regression with adjusted
	# weights (inversely proportional to class frequency)
	lr = LogisticRegression(C = 1.0, penalty = 'l1', tol = 1e-6)
	logReg = lr.fit(train_X, train_Y)
	lr_pred = lr.predict(test_X)
	print(classification_report(test_Y, lr_pred, target_names=["Invalid", "Valid"], digits=4))

	print("Training Random Forest classifier")

	rf = RandomForestClassifier(n_estimators=100)
	scores = cross_val_score(rf, train_X, train_Y, n_jobs=args.threads, scoring="f1", cv=5)
	print("Random Forest F1 score: " + str(scores.mean()) + " (+/- " + str(scores.std() * 2) + ")")
	rf.fit(train_X, train_Y)
	clf1_y_pred = rf.predict(test_X)
	print(classification_report(test_Y, clf1_y_pred, target_names=["Invalid", "Valid"], digits=4))


main()
