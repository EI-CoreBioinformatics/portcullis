#!/usr/bin/env python3

"""
This python script is intended to merge output from several junction files into one by various means
"""

import collections
import sys
from enum import Enum, unique

from junctools.junction import *

__author__ = "Dan Mapleson"
__copyright__ = "Copyright 2016, Portcullis"
__credits__ = ["Dan Mapleson", "Luca Venturini", "David Swarbreck"]
__license__ = "GPLv3"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"

@unique
class Mode(Enum):
	INTERSECTION = 1
	UNION = 2
	CONSENSUS = 3
	SUBTRACT = 4
	DIFFERENCE = 5

	def multifile(self):
		return self.value <= 3


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

def setops(args):

	mode = Mode[args.mode.name.upper()]

	if mode.multifile():

		min_entry = 0

		if mode == Mode.INTERSECTION:
			min_entry = len(args.input)
		elif mode == Mode.UNION:
			min_entry = 1
		elif mode == Mode.CONSENSUS:
			min_entry = args.min_entry
		else:
			raise "Something strange happened"

		if min_entry <= 0:
			raise "Invalid value for min_entry.  Please enter a value of 2 or more."

		merged = collections.defaultdict(list)

		# Check all input files have the same extension
		last_ext = None
		for f in args.input:
			filename, ext = os.path.splitext(f)
			if last_ext != None:
				if last_ext != ext:
					print("Not all input files have the same extension.", out=sys.stderr)
					exit(1)
			else:
				last_ext = ext

		filename, ext = os.path.splitext(args.output)

		if last_ext != ext:
			print("Output extension is not the same as the input.", out=sys.stderr)
			exit(1)

		print("\t".join(["File","distinct", "total"]))

		for f in args.input:
			counter = 0
			found = set()
			with open(f) as fin:
				for line in fin:
					junc = JuncFactory.create_from_ext(last_ext, use_strand=not args.ignore_strand).parse_line(line, fullparse=False)
					if junc:
						counter += 1
						key = junc.key
						found.add(key)
						merged[key].append(line)

			print("\t".join([f, str(len(found)), str(counter)]))

		print("Found a total of", len(merged), "distinct entries.")

		i = 0
		with open(args.output, "wt") as out:

			description = "Merge of multiple junction files.  Min_Entry: {0}. Score_op: {1}".format(min_entry, args.operator.upper())
			header = ExonJunction.create(last_ext).file_header(description=description)
			print(header, file=out)

			for b in sorted(merged):
				if len(merged[b]) >= min_entry:
					juncs = []
					scores = []
					lefts = []
					rights = []
					for line in merged[b]:
						j = ExonJunction.create(last_ext).parse_line(line)
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

def add_options(parser):

	parser.add_argument("-m", "--min_entry", type=int, default="1",
							  help='''Minimum number of files the entry is require to be in.  0 means entry must be
							  present in all files, i.e. true intersection.  1 means a union of all input files''')
	parser.add_argument("--operator", default="total",
							  help="Operator to use for calculating the score in the merged file.  Options: [min, max, sum, mean].  Applicable to intersection, union and consensus modes.")
	parser.add_argument("-o", "--output", required=True, help="Merged output file")
	parser.add_argument("-p", "--prefix", default="junc_merged",
							  help="Prefix to apply to name column in BED output file")
	parser.add_argument("-is", "--ignore_strand", action='store_true', default=False,
								help="Whether or not to ignore strand when creating a key for the junction")
	parser.add_argument("mode", help="Set operation to apply.  Available options: [intersection, union, consensus, subtract, diff]")
	parser.add_argument("input", nargs="+", help="List of junction files to merge (must all be the same type)")