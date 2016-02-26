#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse
import copy
from bed12 import *

def calcOp(op, newval, currval) :
	if op == "max":
		return max(newval, currval)
	elif op == "total":
		return newval + currval
	else:
		# Default to total
		return newval + currval


parser=argparse.ArgumentParser("Script to merge multiple BED files, using target sequence, and thick_start and thick_end values to generate set keys.  (Only entering a single file here essentially copies the file).  Merged BED entries contain the maximum score found across all input files.")
parser.add_argument("input", nargs="+", help="List of BED files to intersect")
parser.add_argument("-m", "--min_entry", type=int, default="1", help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.  1 means a union of all input files")
parser.add_argument("--operator", default="total", help="Operator to use for calculating the score in the merged BAM file.  Options: [max, total]")
parser.add_argument("-o", "--output", required=True, help="Output BED file")
parser.add_argument("-s", "--ignore_strand", default=True, help="Additionally, use strand information for calculating match in bed files")
parser.add_argument("-p", "--prefix", default="Junc_intersected", help="Prefix to apply to name column in BED output file")
args=parser.parse_args()

min_entry = args.min_entry if args.min_entry > 0 else len(args.input)

print("Requested " + str(len(args.input)) + " input BED files to be intersected", file=sys.stderr)
print("Each entry must be found in at least " + str(min_entry) + " files to be in output", file=sys.stderr)

bed_intersect = {}
for bf in args.input:
	bed_set = loadbed(bf, not args.ignore_strand, False)
	print("Input file: " + bf + "; contains " + str(len(bed_set)) + " entries", file=sys.stderr)
	for b in bed_set:
		if b.key not in bed_intersect:
			bed_intersect[b.key] = [b]
		else:
			bed_intersect[b.key].append(b)

print("Found a total of " + str(len(bed_intersect)) + " distinct entries", file=sys.stderr)

bed_out = list()
i = 1
for b in sorted(bed_intersect):
	if len(bed_intersect[b]) >= min_entry:
		bo = copy.deepcopy(bed_intersect[b][0])
		bo.name = args.prefix + "_" + str(i)
		i += 1
		bo.score = 0
		for r in bed_intersect[b]:
			bo.score = calcOp(args.operator, r.score, bo.score)
		bed_out.append(bo)

print("Filtering out " + str(len(bed_intersect) - len(bed_out)) + " entries", file=sys.stderr)
print("Output file contains " + str(len(bed_out)) + " entries", file=sys.stderr)

saveList(args.output, bed_out, "Merge of multiple BED files.  Min_Entry: " + str(min_entry) + ". Score: MAX")

print("Saved file to: " + args.output, file=sys.stderr)
