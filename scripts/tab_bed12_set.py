#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse
from bed12 import *

   
parser=argparse.ArgumentParser("Script to analyse the results of portcullis when run on data from a model organism (i.e. it has a high quality annotated genome.")
parser.add_argument("-a", required=True,
                    help="Portcullis Junction Tab file.  Details from this file will be used in the output.")
parser.add_argument("-b", required=True,
                    help="Bed file")
parser.add_argument("-o", required=True,
                    help="Output tab file")
parser.add_argument("-t", "--tophat", action='store_true', default=False,
                    help="Bed file is from tophat, compensate for maximal overhangs")
parser.add_argument("-s", "--ignore_strand", action='store_true', default=False,
                    help="Use strand information in bed files")
parser.add_argument("-f", "--function", required=True,
                    help="Function to apply: intersect, subtract")
args=parser.parse_args()

mode = -1

if args.function == "intersect":
    mode = 0
elif args.function == "subtract":
    mode = 1

else:
    print("Unrecognised function")
    exit(3)

if mode == -1:
    exit(1)


else:
    # Load second bed file as a set
    bed = loadbed(args.b, not args.ignore_strand, args.tophat)

    filtertab(args.a, args.o, bed, mode, args.ignore_strand)
