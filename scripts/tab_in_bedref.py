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
	parser.add_argument("-f", "--filter", action='store_true', default=False, help="Whether to filter tab file")
	args = parser.parse_args()

	# X should contain a matrix of features derived from the portcullis tab file
	# y should contain the labels (0 not a valid junction, 1 a valid junction).  Confirmed with the reference.

	# Load tab file and produce matrix

	header = ""

	with open(args.input[0]) as f:
		# Skip header
		header = f.readline()
	f.close()

	# b, t = tab.loadtab(i)
	# bed.extend(b)
	# tabs.extend(t)
	# print ("Loaded " + str(len(b)) + " entries from: " + i, file=sys.stderr)
	#    print ("# tab entries: " + str(len(tabs)) + " from " + str(len(args.input)) + " input files", file=sys.stderr)


	# Load reference and add labels
	ref = bed12.loadbed(args.reference, False, False)
	print("# ref entries: " + str(len(ref)), file=sys.stderr)

	res = open(args.output, "w")

	filtin = None
	filtout = None

	if args.filter:
		filtin = open(args.output + ".in.tab", "w")
		filtout = open(args.output + ".out.tab", "w")
		filtin.write(header)
		filtout.write(header)

	nbentries = 0
	for tf in args.input:
		with open(tf) as f:

			# Skip header
			h = f.readline()

			for line in f:

				cleanline = line.strip()

				if not cleanline == "":
					b = bed12.BedEntry.create_from_tabline(line, False, False)
					nbentries += 1
					if b in ref:
						print("1", file=res)
						if args.filter:
							print(line, file=filtin, end="")
					else:
						print("0", file=res)
						if args.filter:
							print(line, file=filtout, end="")

	res.close()

	if args.filter:
		filtin.close()
		filtout.close()

	print("Found ", nbentries, " tab entries in ", len(args.input), " input files")


main()
