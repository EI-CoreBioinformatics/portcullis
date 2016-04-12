#!/usr/bin/env python3

import sys
import os
import argparse
import bed12
import tab


def main():
	parser = argparse.ArgumentParser("Script to produce a list labelling whether each entry in the portcullis tab input belong to the bed reference or not")
	parser.add_argument("input", nargs="+", help="The tab file produce by portcullis")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-o", "--output", default="bedref_out.labels", help="Output prefix for output files")
	args = parser.parse_args()

	# X should contain a matrix of features derived from the portcullis tab file
	# y should contain the labels (0 not a valid junction, 1 a valid junction).  Confirmed with the reference.

	# Load tab file and produce matrix

	header = ""

	with open(args.input[0]) as f:
		# Skip header
		header = f.readline()
	f.close()

	# Load reference and add labels
	ref = bed12.loadbed(args.reference, False, False)
	print("# ref entries: " + str(len(ref)), file=sys.stderr)

	res = open(args.output, "w")

	nbentries = 0
	nb_pos = 0
	nb_neg = 0
	for tf in args.input:
		with open(tf) as f:

			# Skip header
			h = f.readline()

			for line in f:

				cleanline = line.strip()

				if not cleanline == "":
					b = bed12.BedEntry.create_from_line(line, False, False)
					nbentries += 1
					if b in ref:
						print("1", file=res)
						nb_pos += 1

					else:
						print("0", file=res)
						nb_neg += 1


	res.close()

	print("Found", nbentries, "bed entries in", len(args.input), "input files")
	print("Detected", nb_pos, "positive and", nb_neg, "negative entries")


main()
