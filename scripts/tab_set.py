#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse
from tab import *

parser = argparse.ArgumentParser("Script to perform set operations on portcullis tab files.")
parser.add_argument("-a", required=True,
					help="First tab file.  Details from this file will be used in the output.")
parser.add_argument("-b", required=True, help="Second tab file")
parser.add_argument("-o", required=True, help="Output file")
parser.add_argument("-f", "--function", required=True,
					help="Function to apply: intersect, subtract, union")
args = parser.parse_args()

mode = -1

if args.function == "intersect":
	mode = 0
elif args.function == "subtract":
	mode = 1
elif args.function == "union":
	mode = 2
else:
	print("Unrecognised function")
	exit(3)

if mode == -1:
	exit(1)

if mode == 2:

	exit(2)


else:
	# Load second bed file as a set
	tab = loadtabasset(args.b)

	print("Found " + str(len(tab)) + " entries in B file")

	filtertab(args.a, args.o, tab, mode)
