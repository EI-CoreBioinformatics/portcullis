#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse

def makekey(line, usestrand, tophat) :
    words = line.split()
    overhang = words[10]
    overhang_parts = overhang.split(",")
    lo = int(overhang_parts[0])
    ro = int(overhang_parts[1])
    chr = words[0]
    start = str(int(words[6]) + lo) if tophat else words[6]
    end = str(int(words[7]) - ro - 1) if tophat else str(int(words[7]) - 1)
    strand = words[5]
    key = chr + "_" + start + "_" + end
    if usestrand:
        key += "_" + strand
    return key

def makekeyfromtab(line, usestrand) :
    words = line.split()
    chr = words[2]
    start = words[4]
    end = words[5]
    strand = words[11]
    key = chr + "_" + start + "_" + end
    if usestrand:
        key += "_" + strand
    return key
    

def loadbed(filepath, usestrand, tophat) :
    with open(filepath) as f:
        index = 0
        items = set()
        for line in f:
            if index > 0:
                key = makekey(line, usestrand, tophat)
                items.add(key)
            index += 1
    if len(items) != index - 1 :
        print ("non unique items in bed file " + filepath)
    return items



   
parser=argparse.ArgumentParser("Script to analyse the results of portcullis when run on data from a model organism (i.e. it has a high quality annotated genome.")
parser.add_argument("-a", required=True,
                    help="Portcullis Junction Tab file.  Details from this file will be used in the output.")
parser.add_argument("-b", required=True,
                    help="Bed file")
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

    index = 0;
    with open(args.a) as f:

        print(f.readline(), end="")
        for line in f:

            line = line.strip()
            if line != "":
                key = makekeyfromtab(line, not args.ignore_strand)
                if mode == 0:
                    if key in bed:
                        print(line)
                elif mode == 1:
                    if key not in bed:
                        print(line)
