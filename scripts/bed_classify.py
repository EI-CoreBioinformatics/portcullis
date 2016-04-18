#!/usr/bin/env python3

import sys
import os
import argparse
import bed12
from performance import Performance


def main():
	parser = argparse.ArgumentParser("Script to produce a list labelling whether each entry in the portcullis tab input belong to the bed reference or not.  Please note that this will only work with reference files that have distinct entries.")
	parser.add_argument("input", help="The filtered bed file to assess")
	parser.add_argument("-r", "--reference", required=True, help="The input bed file prior to filtering")
	parser.add_argument("-l", "--labels", required=True, help="The labels which should match entries in reference")
	parser.add_argument("-o", "--output", default="bedref_out.labels", help="Output prefix for output files")
	args = parser.parse_args()

	# Load reference and labels, divide into tp and tn
	rp = set()
	rn = set()

	labs = open(args.labels, "r")
	refs = open(args.reference, "r")

	# read header line from bed file
	header = refs.readline()

	# Now we should have the same number of lines in both files
	line = 1
	while 1:
		ref_line = refs.readline()
		lab_line = labs.readline()
		if not ref_line and not lab_line: break
		if (lab_line and not ref_line) or (not lab_line and ref_line): print("ERROR: reference file and labels file have a different number of entries.", file=sys.stderr); exit(1)
		ref_line = ref_line.strip()
		lab_line = lab_line.strip()
		if lab_line == "1":
			rp.add(bed12.BedEntry.create_from_line(ref_line, False, False))
		elif lab_line == "0":
			rn.add(bed12.BedEntry.create_from_line(ref_line, False, False))
		else:
			print("ERROR: Label file contains an entry that is not either \"0\" or \"1\" at line:", line, file=sys.stderr); exit(1)
		line += 1
	labs.close()
	refs.close()

	print("Reference contains", len(rp), "positive and", len(rn), "negative entries.")

	p = set()
	with open(args.input) as f:

		# Skip header
		h = f.readline()

		for line in f:

			cleanline = line.strip()

			if not cleanline == "":
				p.add(bed12.BedEntry.create_from_line(line, False, False))

	print("Input contains", len(p), "entries")

	tp = p & rp
	fp = p & rn
	fn = rp - p
	tn = rn - p

	perf = Performance(tp = len(tp), fp = len(fp), fn = len(fn), tn = len(tn))
	print(Performance.longHeader())
	print(perf.longStr())

main()
