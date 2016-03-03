#!/usr/bin/env python3

import os
import argparse
import bed12
from sklearn.cross_validation import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import classification_report

import numpy as np
import matplotlib.pyplot as plt

class TabEntry:
	chrom = ""
	start = 0
	end = 0
	left = 0
	right = 0
	strand = "?"
	#M1 = "N"
	M2 = 0
	M3 = 0
	M4 = 0
	#M5 = 0
	#M6 = 0
	#M7 = 0
	M8 = 0
	M9 = 0
	M10 = 0
	M11 = 0.0
	M12 = 0
	M13 = 0
	M14 = 0
	#M15 = 0.0

	def __init__(self):
		self.data = []

	def __str__(self):

		line = [ self.chrom, self.start, self.end, self.left, self.right, self.strand,
			 self.M2, self.M3, self.M4, self.M8, self.M9, self.M10, self.M11, self.M12, self.M13, self.M14
			 ]
		return "\t".join([str(_) for _ in line])

	def makeMatrixRow(self):
		return [ self.M2, self.M3, self.M4, self.M8, self.M9, self.M10, self.M11, self.M12, self.M13, self.M14 ]

	@staticmethod
	def features():
		return [ "M2-nb-reads", "M3-nb_dist_aln", "M4-nb_rel_aln", "M8-max_min_anc", "M9-dif_anc", "M10-dist_anc",
				 "M11-entropy", "M12-maxmmes", "M13-hammping5p", "M14-hamming3p"]
	@staticmethod
	def featureAt(index):
		return TabEntry.features()[index]

	@staticmethod
	def sortedFeatures(indicies):
		f = TabEntry.features()
		s=[]
		for i in indicies:
			s.append(f[i])
		return s

	@staticmethod
	def nbMetrics():
		return 10

	@staticmethod
	def create_from_tabline(line):

		b = TabEntry()

		parts = line.split("\t")

		b.chrom = parts[2]
		b.left = int(parts[6])
		b.right = int(parts[7])
		b.strand = parts[12]
		b.start = int(parts[4])
		b.end = int(parts[5])

		b.M2 = int(parts[14])
		b.M3 = int(parts[15])
		b.M4 = int(parts[16])
		b.M8 = int(parts[20])
		b.M9 = int(parts[21])
		b.M10 = int(parts[22])
		b.M11 = float(parts[23])
		b.M12 = int(parts[24])
		b.M13 = int(parts[25])
		b.M14 = int(parts[26])

		return b



def loadtab(tabfile):

	bed=list()
	tab=list()
	with open(tabfile) as f:
		# Skip header
		f.readline()

		for line in f:
			line.strip()
			if len(line) > 1:
				bed.append(bed12.BedEntry.create_from_tabline(line, False, False))
				tab.append(TabEntry.create_from_tabline(line))


	return bed, tab



def main():

	parser = argparse.ArgumentParser("Script to build a random forest decision tree")
	parser.add_argument("input", nargs="+", help="The tab file produce by portcullis")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-t", "--threads", type=int, default="1", help="The number of threads to use")
	parser.add_argument("--test", help="Test the classifier against this file")
	parser.add_argument("-o", "--output", required=True, help="The output prefix")
	args = parser.parse_args()

	# X should contain a matrix of features derived from the portcullis tab file
	# y should contain the labels (0 not a valid junction, 1 a valid junction).  Confirmed with the reference.

	# Load tab file and produce matrix
	bed=[]
	tab=[]
	for i in args.input:
		b, t = tab.loadtab(i)
		bed.extend(b)
		tab.extend(t)
		print ("Loaded " + str(len(b)) + " entries from: " + i)
	print ("# tab entries: " + str(len(tab)) + " from " + str(len(args.input)) + " input files")


	# Load reference and add labels
	ref = bed12.loadbed(args.reference, False, False)
	print ("# ref entries: " + str(len(ref)))

	in_juncs = 0
	out_juncs = 0
	X=np.zeros( (len(tab), TabEntry.nbMetrics()) )
	y=list()
	for i in range(0, len(bed)):
		b = bed[i]
		X[i] = tab[i].makeMatrixRow()
		if b in ref:
			in_juncs += 1
			y.append(1)
		else:
			out_juncs += 1
			y.append(0)

	print ("In:" + str(in_juncs))
	print ("Out:" + str(out_juncs))

	# Load test data
	test_b, test_t = loadtab(args.test)
	test_X=np.zeros( (len(test_t), TabEntry.nbMetrics()) )
	test_y=[]
	for i in range(0, len(test_t)):
		b = test_b[i]
		test_X[i] = test_t[i].makeMatrixRow()
		if b in ref:
			test_y.append(1)
		else:
			test_y.append(0)


	print("Training Random Forest classifier")

	clf1 = RandomForestClassifier(n_estimators=40)
	scores = cross_val_score(clf1, X, y, n_jobs=args.threads, scoring="f1", cv=10)
	print("Random Forest F1 score: " + str(scores.mean()) + " (+/- " + str(scores.std() * 2))
	clf1.fit(X, y)
	clf1_y_pred = clf1.predict(test_X)
	print(classification_report(test_y, clf1_y_pred, target_names=["Invalid", "Valid"]))


	#print("Training SVM (with RBF) classifier")
	#clf2 = SVC()
	#scores = cross_val_score(clf2, X, y, n_jobs=args.threads, scoring="f1")
	#print("SVM Mean score: " + str(scores.mean()))
	#clf2.fit(X, y)
	#clf2_y_pred = clf2.predict(test_X)
	#print(classification_report(test_y, clf2_y_pred, target_names=["0", "1"]))


	importances = clf1.feature_importances_
	std = np.std([tree.feature_importances_ for tree in clf1.estimators_], axis=0)
	indices = np.argsort(importances)[::-1]




	# Print the feature ranking
	print("Feature ranking:")

	for f in range(X.shape[1]):
		print("%d. feature %s (%f)" % (f + 1, TabEntry.featureAt(indices[f]), importances[indices[f]]))

	# Plot the feature importances of the forest
	plt.figure()
	plt.title("Feature importances")
	plt.bar(range(X.shape[1]), importances[indices],
		   color="r", yerr=std[indices], align="center")
	plt.xticks(range(X.shape[1]), TabEntry.sortedFeatures(indices))
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	plt.xlim([-1, X.shape[1]])
	plt.tight_layout()

	plt.savefig(args.output + ".png")


	# Print

	
main()