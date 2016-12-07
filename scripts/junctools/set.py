#!/usr/bin/env python3

"""
This python script is intended to merge output from several junction files into one by various means
"""

import argparse
import collections

from .junction import *

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
	SYMMETRIC_DIFFERENCE = 5
	FILTER = 6
	IS_SUBSET = 7
	IS_SUPERSET = 8
	IS_DISJOINT = 9

	def multifile(self):
		return self.value <= 3

	def makes_output(self):
		return self.value <= 6

	def is_test(self):
		return self.value >= 7

	def needs_consistent_ext(self):
		return self.multifile() or self.name == Mode.SYMMETRIC_DIFFERENCE.name


@unique
class CalcOp(Enum):
	MIN = 1
	MAX = 2
	SUM = 3
	MEAN = 4

	def execute(self, vals):
		if self.value == CalcOp.MAX.value:
			return max(vals)
		elif self.value == CalcOp.MIN.value:
			return min(vals)
		elif self.value == CalcOp.SUM.value:
			return sum(vals)
		elif self.value == CalcOp.MEAN.value:
			return sum(vals) / float(len(vals))
		else:
			raise ValueError("CalcOp Error - Should never happen")



def setops(args):
	mode = Mode[args.mode.upper()]

	min_entry = 1

	if len(args.input) < 2:
		raise ValueError("We require a least two input files")

	if not mode.multifile() and len(args.input) > 2:
		raise ValueError("We don't support more than two input files for this mode")

	if mode.multifile():
		if mode == Mode.INTERSECTION:
			min_entry = len(args.input)
		elif mode == Mode.UNION:
			min_entry = 1
		elif mode == Mode.CONSENSUS:
			min_entry = args.min_entry
		else:
			raise ValueError("Something strange happened")

	if min_entry <= 0:
		raise ValueError("Invalid value for min_entry.  Please enter a value of 2 or more.")

	# Check all input files have the same extension if required
	last_ext = None
	if mode.needs_consistent_ext():
		for f in args.input:
			filename, ext = os.path.splitext(f)
			if last_ext != None:
				if last_ext != ext:
					raise ValueError("Not all input files have the same extension.")
			else:
				last_ext = ext
	else:
		filename, ext = os.path.splitext(args.input[0])
		last_ext = ext

	if mode.makes_output():

		if not args.output:
			raise ValueError("Set operation requested that produces output.  Please set the '--output' option.")

		filename, ext = os.path.splitext(args.output)

		if last_ext != ext:
			raise ValueError("Output extension is not the same as the input.")

	print("junctools set", mode.name.lower())

	if mode.multifile():
		merged = collections.defaultdict(list)

		print("\t".join(["File", "distinct", "total"]))

		for f in args.input:
			counter = 0
			found = set()
			with open(f) as fin:
				for line in fin:
					junc = JuncFactory.create_from_ext(last_ext, use_strand=not args.ignore_strand).parse_line(line,
																											   fullparse=False)
					if junc:
						counter += 1
						key = junc.key
						if not key in found:
							found.add(key)
						merged[key].append(line)

			print("\t".join([f, str(len(found)), str(counter)]))

		print()
		print("The union of all files contains", len(merged), "distinct entries.")

		i = 0
		with open(args.output, "wt") as out:

			description = "Set operation on junction files. Mode: {0};  Min_Entry: {1}; Score_op: {2}".format(mode.name,
																											  min_entry,
																											  args.operator.upper())
			header = JuncFactory.create_from_ext(last_ext).file_header(description=description)
			print(header, file=out)

			calcop = CalcOp[args.operator.upper()]

			for b in sorted(merged):
				nb_samples = len(merged[b])
				if nb_samples >= min_entry:
					juncs = []
					scores = []
					lefts = []
					rights = []
					for line in merged[b]:
						j = JuncFactory.create_from_ext(last_ext).parse_line(line)
						juncs.append(j)
						scores.append(j.score)
						lefts.append(j.left)
						rights.append(j.right)

					merged_junc = juncs[0]
					merged_junc.name = "{prefix}_{i}".format(prefix=args.prefix, i=i)
					merged_junc.score = calcop.execute(scores)
					merged_junc.left = CalcOp.MIN.execute(lefts)
					merged_junc.right = CalcOp.MAX.execute(rights)

					if type(merged_junc) is TabJunction:
						merged_junc.setNbSamples(nb_samples)

					i += 1

					print(merged_junc, file=out)

			if mode != Mode.UNION:
				print("Filtered out", len(merged) - i, "entries")
			print("Output file", args.output, "contains", i, "entries.")

	else:
		if mode.makes_output():


			with open(args.output, "wt") as out:
				description = "Set operation on junction files. Mode: {0}".format(mode.name)
				header = JuncFactory.create_from_file(args.input[0]).file_header(description=description)
				print(header, file=out)

				out_count = 0
				if mode == Mode.SUBTRACT or mode == Mode.FILTER:
					print("Loading second input file into a set")
					ref, entries = Junction.createJuncSet(args.input[1], use_strand=not args.ignore_strand, fullparse=False)
					print("\t".join(["File", "Total", "Distinct"]))
					print("\t".join([args.input[1], str(entries), str(len(ref))]))

					with open(args.input[0]) as f:
						for line in f:
							junc = JuncFactory.create_from_file(args.input[0], use_strand=not args.ignore_strand).parse_line(line,
																											fullparse=False)
							if junc and ((mode == Mode.SUBTRACT and not junc.key in ref) or (mode == Mode.FILTER and junc.key in ref)):
								print(line.rstrip(), file=out)
								out_count += 1

				elif mode == Mode.SYMMETRIC_DIFFERENCE:
					print("Loading input files into sets")
					print("\t".join(["File", "Total", "Distinct"]))
					sets = []
					for f in args.input:
						juncs, entries = Junction.createJuncSet(f, use_strand=not args.ignore_strand, fullparse=False)
						sets.append(juncs)
						print("\t".join([f, str(entries), str(len(juncs))]))
					print()

					with open(args.input[0]) as f:
						for line in f:
							junc = JuncFactory.create_from_file(args.input[0], use_strand=not args.ignore_strand).parse_line(line,
																											fullparse=False)
							if junc:
								if not junc.key in sets[1]:
									print(line.rstrip(), file=out)
									out_count += 1

					with open(args.input[1]) as f:
						for line in f:
							junc = JuncFactory.create_from_file(args.input[1], use_strand=not args.ignore_strand).parse_line(line,
																											fullparse=False)
							if junc:
								if not junc.key in sets[0]:
									print(line.rstrip(), file=out)
									out_count += 1


				print()
				print("Output contains ", out_count, " junctions")
				print("Output saved to", args.output)

		elif mode.is_test():

			print()
			print("\t".join(["File", "Total", "Distinct"]))
			sets = []
			for f in args.input:
				juncs, entries = Junction.createJuncSet(f, use_strand=not args.ignore_strand, fullparse=False)
				sets.append(juncs)
				print("\t".join([f, str(entries), str(len(juncs))]))
			print()

			res = None

			if mode == Mode.IS_SUBSET:
				res = sets[0].issubset(sets[1])
			elif mode == Mode.IS_SUPERSET:
				res = sets[0].issuperset(sets[1])
			elif mode == Mode.IS_DISJOINT:
				res = sets[0].isdisjoint(sets[1])

			print("True" if res else "False")

		else:
			raise "Unknown mode"


def add_options(parser):
	parser.formatter_class = argparse.RawTextHelpFormatter

	parser.add_argument("-m", "--min_entry", type=int, default="1",
						help='''Minimum number of files the entry is require to be in.  0 means entry must be
present in all files, i.e. true intersection.  1 means a union of all input files''')
	parser.add_argument("--operator", default="sum",
						help='''Operator to use for calculating the score in the merged file.
This option is only applicable to 'intersection', 'union' and 'consensus' modes.
Available values:
 - min
 - max
 - sum
 - mean''')
	parser.add_argument("-o", "--output",
						help="Output junction file.  Required for operations that produce an output file.")
	parser.add_argument("-p", "--prefix", default="junc_merged",
						help="Prefix to apply to name column in BED output file")
	parser.add_argument("-is", "--ignore_strand", action='store_true', default=False,
						help="Whether or not to ignore strand when creating a key for the junction")
	parser.add_argument("mode", help='''Set operation to apply.  See above for details.  Available options:
 - intersection
 - union
 - consensus
 - subtract
 - filter
 - symmetric_difference
 - is_subset
 - is_superset
 - is_disjoint
''')
	parser.add_argument("input", nargs="+", help="List of junction files to merge (must all be the same file format)")
