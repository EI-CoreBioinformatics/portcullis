#!/usr/bin/env python3

"""
This python script is intended to merge output from several junction files into one by various means
"""

import sys
import os
import argparse
import collections
from portcullis_junction import *

__author__ = "Dan Mapleson"
__copyright__ = "Copyright 2016, Portcullis"
__credits__ = ["Dan Mapleson", "Luca Venturini", "David Swarbreck"]
__license__ = "GPLv3"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"

def calc_op(op, vals):
	if op == "max":
		return max(vals)
	elif op == "min":
		return min(vals)
	elif op == "sum":
		return sum(vals)
	elif op == "mean":
		return sum(vals) / float(len(vals))
	else:
		raise ValueError

def main():

	parser = argparse.ArgumentParser("This script can merge output from several portcullis runs into one.\
		This supports BED12 and TAB formats that contain exon anchors.\
		Please note for TAB format that while anchors will be extended the metrics used to a merged junction will reflect\
		only the junction from the first sample and are therefore not recalculated based on the merged information")

	parser.add_argument("-m", "--min_entry", type=int, default="1",
							help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.  1 means a union of all input files")
	parser.add_argument("--operator", default="total",
							help="Operator to use for calculating the score in the merged file.  Options: [min, max, sum, mean]")
	parser.add_argument("-o", "--output", required=True, help="Merged output file")
	parser.add_argument("-p", "--prefix", default="junc_merged",
							help="Prefix to apply to name column in BED output file")
	parser.add_argument("input", nargs="+", help="List of BED or TAB files to merge (must all be the same type)")

	args = parser.parse_args()

	min_entry = args.min_entry if args.min_entry > 0 else len(args.input)

	merged = collections.defaultdict(list)

	# Check all input files have the same extension
	last_ext = None
	for f in args.input:
		filename, ext = os.path.splitext(f)
		if not last_ext == None:
			if not last_ext == ext:
				print("Not all input files have the same extension.", out=sys.stderr)
				exit(1)
		else:
			last_ext = ext


	for f in args.input:
		counter = 0
		found = set()
		with open(f) as fin:
			for line in fin:
				junc = create_exon_junction(last_ext)
				res = junc.parse_line(line, fullparse=False)
				if res is None:
					# Skipping header
					continue
				counter += 1
				key = junc.key
				found.add(key)
				merged[key].append(line)

		if len(found) < counter:
			print("WARNING: File", f, " had duplicated entries (", counter, "total entries,", len(found), "distinct entries)", out=sys.stderr)
		print("Input file", f, "contains", len(found), "entries.")

	print("Found a total of", len(merged), "distinct entries.")

	i = 0
	with open(args.output, "wt") as out:

		description = "Merge of multiple junction files.  Min_Entry: {0}. Score_op: {1}".format(min_entry, args.operator.upper())
		header = create_exon_junction(last_ext).file_header(description=description)

		print(header, file=out)
		for b in sorted(merged):
			if len(merged[b]) >= min_entry:
				juncs = []
				scores = []
				lefts = []
				rights = []
				for line in merged[b]:
					j = create_exon_junction(last_ext).parse_line(line)
					juncs.append(j)
					scores.append(j.score)
					lefts.append(j.left)
					rights.append(j.right)

				merged_junc = juncs[0]
				merged_junc.name = "{prefix}_{i}".format(prefix=args.prefix, i=i)
				merged_junc.score = calc_op(args.operator, scores)
				merged_junc.left = calc_op("min", lefts)
				merged_junc.right = calc_op("max", rights)

				i += 1

				print(merged_junc, file=out)

	print("Filtered out", len(merged)-i, "entries")
	print("Output file", args.output, "contains", i, "entries.")

if __name__ == '__main__':
	main()
