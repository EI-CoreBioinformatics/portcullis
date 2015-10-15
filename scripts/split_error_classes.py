#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse
from bed12 import *

parser=argparse.ArgumentParser("Script to analyse the results of portcullis when run on data from a model organism (i.e. it has a high quality annotated genome.")
parser.add_argument("-r", required=True,
                    help="Reference BED file")
parser.add_argument("-p", required=True,
                    help="BED file containing junctions passing a filter")
parser.add_argument("-f", required=True,
                    help="BED file containing junctions failing a filter")
parser.add_argument("-u", required=True,
                    help="Unfiltered tab file")
parser.add_argument("-o", required=True,
                    help="Output prefix for output files")
parser.add_argument("-s", "--ignore_strand", action='store_true', default=False,
                    help="Use strand information in bed files")
args=parser.parse_args()


usesstrand = not args.ignore_strand

# Load reference bed file as a set
bedRef = loadbed(args.r, usesstrand, False)

# Filenames
tpfile = args.o + ".TP.bed"
fnfile = args.o + ".FN.bed"
fpfile = args.o + ".FP.bed"
tnfile = args.o + ".TN.bed"

# Intersect reference with passed data - TP
filterbed(args.p, bedRef, 0, usesstrand, False, tpfile)

# Intersect reference with failed data - FN
filterbed(args.f, bedRef, 0, usesstrand, False, fnfile)

# Subtract reference from passed data - FP
filterbed(args.p, bedRef, 1, usesstrand, False, fpfile)

# Subtract reference from failed data - TN
filterbed(args.f, bedRef, 1, usesstrand, False, tnfile)


bedtp = loadbed(tpfile, usesstrand, False)
bedfn = loadbed(fnfile, usesstrand, False)
bedfp = loadbed(fpfile, usesstrand, False)
bedtn = loadbed(tnfile, usesstrand, False)

filtertab(args.u, args.o + ".TP.tab", bedtp, 0, usesstrand)
filtertab(args.u, args.o + ".FN.tab", bedfn, 0, usesstrand)
filtertab(args.u, args.o + ".FP.tab", bedfp, 0, usesstrand)
filtertab(args.u, args.o + ".TN.tab", bedtn, 0, usesstrand)
