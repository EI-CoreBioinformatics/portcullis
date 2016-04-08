#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import numpy as np
from sklearn import cross_validation
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import scale
from mpl_toolkits.mplot3d import Axes3D

import bed12
import tab


def main():
	parser = argparse.ArgumentParser("Script to build a random forest decision tree")
	parser.add_argument("input", nargs="+", help="The tab file produce by portcullis")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	args = parser.parse_args()

	inputstr = ", ".join(args.input)
	# X should contain a matrix of features derived from the portcullis tab file
	# y should contain the labels (0 not a valid junction, 1 a valid junction).  Confirmed with the reference.

	# Load tab file and produce matrix
	tabdata = []
	bed = []
	for i in args.input:
		b, t = tab.loadtab(i)
		tabdata.extend(t)
		bed.extend(b)
		print("Loaded " + str(len(b)) + " entries from: " + i)
	print("# tab entries: " + str(len(tabdata)) + " from " + str(len(args.input)) + " input files")

	# Load reference and add labels
	ref = bed12.loadbed(args.reference, False, False)
	print("# ref entries: " + str(len(ref)))

	in_juncs = 0
	out_juncs = 0
	X = np.zeros((len(tabdata), tab.TabEntry.nbMetrics()))
	y = np.zeros((len(tabdata), 1))
	for i in range(0, len(tabdata)):
		b = bed[i]
		X[i] = tabdata[i].makeMatrixRow()
		if b in ref:
			in_juncs += 1
			y[i, 0] = 1
		else:
			out_juncs += 1
			y[i, 0] = 0

	y = y[:, 0]

	print("In:" + str(in_juncs))
	print("Out:" + str(out_juncs))

	print("Running PCA")

	pca = PCA(n_components=10)
	X_r = pca.fit_transform(scale(X))
	print(pca.explained_variance_ratio_)

	fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(27, 8))
	plt.tight_layout(pad=4.0, w_pad=4.0, h_pad=3.0)

	for c, i, target_name in zip("gr", [1, 0], ["genuine", "invalid"]):
		ax1.scatter(X_r[y == i, 0], X_r[y == i, 1], c=c, label=target_name)
	ax1.legend()
	ax1.set_title("PCA")
	ax1.set_xlabel("PC1")
	ax1.set_ylabel("PC2")

	for c, i, target_name in zip("gr", [1, 0], ["genuine", "invalid"]):
		ax2.scatter(X_r[y == i, 0], X_r[y == i, 1], c=c, label=target_name)
	ax2.legend()
	ax2.set_title("Zoomed PCA")
	ax2.set_xlabel("PC1")
	ax2.set_ylabel("PC2")
	ax2.set_xlim(-5, 10)
	ax2.set_ylim(-2, 10)


	n = len(X_r)
	kf_5 = cross_validation.KFold(n, n_folds=5, shuffle=True, random_state=2)

	regr = LogisticRegression()
	acc = []

	score = cross_validation.cross_val_score(regr, np.ones((n, 1)), y.ravel(), cv=kf_5,
												  scoring='f1').mean()
	acc.append(score)

	for i in np.arange(1, 11):
		score = cross_validation.cross_val_score(regr, X_r[:, :i], y.ravel(), cv=kf_5,
													  scoring='f1').mean()
		acc.append(score)

	print(acc)

	ax3.plot([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], acc[0:11], '-v')
	ax3.set_title('Logistic regression using 5-fold CV')
	ax3.set_xlabel('Number of principal components in regression')
	ax3.set_ylabel('Accuracy')

	ax4.plot([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], acc[1:11], '-v')
	ax4.set_title('Logistic regression excluding initial')
	ax4.set_xlabel('Number of principal components in regression')
	ax4.set_ylabel('Accuracy')
	ax4.set_xlim(0, 10)

	plt.show()

main()
