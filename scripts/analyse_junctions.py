#!/usr/bin/env python3

import os
import argparse
import bed12


class PEntry:
	junc = ""
	aligner = ""
	tp = 0
	tn = 0
	fp = 0
	fn = 0
	sum = 0
	sen = 0.0
	spc = 0.0
	prc = 0.0
	acc = 0.0
	f1 = 0.0

	def __init__(self):
		self.data = []

	def __str__(self):
		return self.junc + "\t" + self.aligner \
			   + "\t" + str(self.tp) \
			   + "\t" + str(self.fp) \
			   + "\t" + str(self.fn) \
			   + "\t" + format(self.sen, '.2f') + "\t" + format(self.prc, '.2f') \
			   + "\t" + format(self.f1, '.2f')

	def calc_scores(self):
		self.sum = self.tp + self.fp + self.tn + self.fn
		self.sen = (float(self.tp) / float(self.tp + self.fn)) * 100.0 if self.tp + self.fn > 0 else 0.0
		self.spc = (float(self.tn) / float(self.tn + self.fp)) * 100.0 if self.tn + self.fp > 0 else 0.0
		self.prc = (float(self.tp) / float(self.tp + self.fp)) * 100.0 if self.tp + self.fp > 0 else 0.0
		self.acc = (
					   float(self.tp + self.tn) / float(
						   self.tp + self.fp + self.tn + self.fn)) * 100.0 if self.sum > 0 else 0.0
		self.f1 = (float(2 * self.tp) / float(
			2 * self.tp + self.fp + self.fn)) * 100.0 if self.tp + self.fp + self.fn > 0 else 0.0

	@staticmethod
	def header():
		return "JuncTool\tAligner\tTP\tFP\tFN\tREC\tPRC\tF1"


'''
    def __str__(self):
        return self.junc + "\t" + self.aligner \
                        + "\t" + str(self.tp) \
                        + "\t" + str(self.tn) \
                        + "\t" + str(self.fp) \
                        + "\t" + str(self.fn) \
                        + "\t" + format(self.sen, '.2f') + "\t" + format(self.spc, '.2f') + "\t" + format(self.prc, '.2f') \
                        + "\t" + format(self.acc, '.2f') + "\t" + format(self.f1, '.2f')
'''


def main():
	parser = argparse.ArgumentParser("Script to create the Venn Plots from BED files")
	parser.add_argument("input", nargs='+', help="The directory containing BED files from pipeline")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-o", "--output", required=True, help="The output prefix")
	args = parser.parse_args()

	ref_bed = bed12.loadbed(args.reference, False, False)
	print("Loaded Reference BED file.  # junctions: " + str(len(ref_bed)))

	# Load all bed files
	bed_data = {}
	aligners = set()
	reads = set()
	junc_analysers = set()
	for bed_file in args.input:
		bed_base = os.path.splitext(os.path.basename(bed_file))[0]
		parts = bed_base.split('-')
		if (not parts[0] == "trinity"):
			aligners.add(parts[0])
			reads.add(parts[1])
			junc_analysers.add(parts[2])
			bed_data[parts[2] + "-" + parts[0]] = bed12.loadbed(bed_file, False, False)
			print("Loaded: " + bed_file + "; # junctions: " + str(len(bed_data[parts[2] + "-" + parts[0]])))

	print("Found these aligners: " + ', '.join(aligners))
	print("Found these reads: " + ', '.join(reads))
	print("Found these junction analysis tools: " + ', '.join(junc_analysers))

	# Build table
	tab = list()
	for a in aligners:
		for j in junc_analysers:
			p = PEntry()
			p.aligner = a
			p.junc = j
			p.tp = len(bed_data[j + "-" + a] & ref_bed)
			p.fp = len(bed_data[j + "-" + a] - ref_bed)
			p.fn = len(ref_bed - bed_data[j + "-" + a])
			p.calc_scores()
			tab.append(p)

	# Output table to disk
	with open(args.output + "-junc_analysis.tab", "w") as tab_out:
		print(PEntry.header(), file=tab_out)
		for p in tab:
			print(p, file=tab_out)


main()
