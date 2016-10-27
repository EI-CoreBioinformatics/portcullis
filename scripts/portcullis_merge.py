#!/usr/bin/env python3

"""
This python script is intended to merge output from several junction files into one by various means
"""

import sys
import argparse
import collections
import logging
import bed12
import tab

__author__ = "Dan Mapleson"
__copyright__ = "Copyright 2016, Portcullis"
__credits__ = ["Dan Mapleson", "Luca Venturini", "David Swarbreck"]
__license__ = "GPLv3"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"


def bed_parse(args, line):
	return bed12.BedEntry.create_from_line(line, False, False)


def tab_parse(args, line):
	return tab.TabEntry.create_from_tabline(line)


def main():
	call_args = sys.argv[1:]

	parser = argparse.ArgumentParser(
		"""This script can merge output from several portcullis runs into one.""")
	subparsers = parser.add_subparsers(
		title="Output type")

	bed_parser = subparsers.add_parser("bed",
									   help="Merges multiple portcullis BED files into one.")
	bed_parser.add_argument("-m", "--min_entry", type=int, default="1",
							help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.  1 means a union of all input files")
	bed_parser.add_argument("--operator", default="total",
							help="Operator to use for calculating the score in the merged BED file.  Options: [max, total]")
	bed_parser.add_argument("-o", "--output", required=True, help="Output BED file")
	bed_parser.add_argument("-l", "--log", default=None, type=argparse.FileType("w"))
	bed_parser.add_argument("-p", "--prefix", default="junc_merged",
							help="Prefix to apply to name column in BED output file")
	bed_parser.add_argument("input", nargs="+", help="List of BED files to merge")

	bed_parser.set_defaults(func=bed_parse)


	tab_parser = subparsers.add_parser("tab",
									   help="Merges multiple portcullis tab files into one.  Will add additional columns to the TAB file describing the results of the merging process")
	tab_parser.add_argument("-m", "--min_entry", type=int, default="1",
							help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.  1 means a union of all input files")
	tab_parser.add_argument("-o", "--output", required=True, help="Output TAB file")
	tab_parser.add_argument("-l", "--log", default=None, type=argparse.FileType("w"))
	tab_parser.add_argument("-s", "--ignore_strand", action='store_true', default=False,
							help="Ignore strand information for calculating match in TAB files")
	tab_parser.add_argument("-p", "--prefix", default="junc_merged",
							help="Prefix to apply to name column in tab output file")
	tab_parser.add_argument("input", nargs="+", help="List of TAB files to merge")
	tab_parser.set_defaults(func=tab_parse)

	args = parser.parse_args(call_args)
	if hasattr(args, "func"):

		logger = logging.getLogger("merger")
		formatter = logging.Formatter(
			"{asctime} - {name} - {filename}:{lineno} - {levelname} - {funcName} - {message}",
			style="{"
		)
		logger.setLevel(logging.INFO)
		if args.log is not None:
			args.log.close()
			args.log = args.log.name
			handler = logging.FileHandler(args.log)
		else:
			handler = logging.StreamHandler()

		handler.setLevel(logging.INFO)
		handler.setFormatter(formatter)
		logger.addHandler(handler)

		min_entry = args.min_entry if args.min_entry > 0 else len(args.input)

		merged = collections.defaultdict(list)

		for f in args.input:
			counter = 0
			found = set()
			with open(f) as fin:
				counter = 0
				found = set()
				for line in fin:
					l = line.strip()
					if not l == "":
						entry = args.func(args, line)
						if entry is None:
							# Skipping header
							continue
						counter += 1
						found.add(entry.key())
						merged[entry.key()].append(line)
						continue

			if len(found) < counter:
				logger.warning("File %s had duplicated entries (%d lines, %d unique entries)",
							   f, counter, len(found))
			logger.info("Input file %s contains %d entries.", f, len(found))

			del found

		logger.info("Found a total of %d distinct entries.", len(merged))

		i = 0
		with open(args.output, "wt") as out:
			description = "Merge of multiple BED files.  Min_Entry: {0}. Score_op: {1}".format(min_entry,
																							   args.operator.upper())
			print("track name=\"junctions\" description=\"{0}\"".format(description), file=out)
			for b in sorted(bed_merge):
				if len(bed_merge[b]) >= min_entry:
					bo = bed12.BedEntry.create_from_line(bed_merge[b][0], not args.ignore_strand, False)
					bo.name = "{prefix}_{i}".format(prefix=args.prefix, i=i)
					i += 1
					for r in bed_merge[b][1:]:
						bo.score = calcOp(args.operator,
										  BedEntry.create_from_line(r, not args.ignore_strand, False).score,
										  bo.score)
					# for num in (8,10,11):
					#     line[num] = ",".join([str(_) for _ in line[num]])
					print(bo, file=out)

		logger.info("Filtered out %d entries", len(merged) - i)
		logger.info("Output file %s contains %d entries.", args.output, i)
		if args.log is not None:
			handler.close()


	else:
		parser.print_help()

if __name__ == '__main__':
	main()
