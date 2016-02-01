#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse
from bed12 import *


def loadbed(filepath, usestrand, tophat) :

    index = 0
    items = set()

    with open(filepath) as f:
        # Skip header
        f.readline()
        for line in f:
            key = makekey(line, usestrand, tophat)
            items.add(key)
            index += 1
    if len(items) != index :
        print ("duplicated items in bed file " + filepath)
    return items


parser=argparse.ArgumentParser("Script to generate the intersection of multiple BED files, using target sequence, and thick_start and thick_end values to generate set keys.  (Only entering a single file here essentially copies the file)")
parser.add_argument("input", nargs="+", help="List of BED files to intersect")
parser.add_argument("-m", "--min_entry", type=int, default="0", help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.")
parser.add_argument("-o", "--output", required=True, help="Output BED file")
parser.add_argument("-s", "--ignore_strand", default=True, help="Additionally, use strand information for calculating match in bed files")
parser.add_argument("-p", "--prefix", default="Junc_intersected", help="Prefix to apply to name column in BED output file")
args=parser.parse_args()

min_entry = args.min_entry if args.min_entry > 0 else len(args.input)

print("Requested " + str(len(args.input)) + " input BED files to be intersected")
print("Each entry must be found in at least " + str(min_entry) + " files to be in output")

bed_intersect = {}
for bf in args.input:
	bed_set = loadbed(bf, not args.ignore_strand, False)
	print("Input file: " + bf + "; contains " + str(len(bed_set)) + " entries")
	for b in bed_set:
		bed_intersect[b] = 1 if not b in bed_intersect else bed_intersect[b] + 1

print("Found a total of " + str(len(bed_intersect)) + " distinct entries")

bed_out = list()
i = 1
for b in bed_intersect:
	if bed_intersect[b] >= min_entry:
		be = BedEntry.create_from_key(b)
		be.name = args.prefix + "_" + str(i)
		i += 1
		bed_out.append(be)

print("Filtering out " + str(len(bed_intersect) - len(bed_out)) + " entries")
print("Output file contains " + str(len(bed_out)) + " entries")

saveList(args.output, bed_out)

print("Saved file to: " + args.output)
