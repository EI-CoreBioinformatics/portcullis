#!/usr/bin/env python3

import sys
import os
import argparse
import bed12
import tab



def main():

    parser = argparse.ArgumentParser("Script to build a random forest decision tree")
    parser.add_argument("input", nargs="+", help="The tab file produce by portcullis")
    parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
    args = parser.parse_args()

    # X should contain a matrix of features derived from the portcullis tab file
    # y should contain the labels (0 not a valid junction, 1 a valid junction).  Confirmed with the reference.

    # Load tab file and produce matrix
    bed=[]
    tabs=[]

    for i in args.input:
        b, t = tab.loadtab(i)
        bed.extend(b)
        tabs.extend(t)
        print ("Loaded " + str(len(b)) + " entries from: " + i, file=sys.stderr)
    print ("# tab entries: " + str(len(tabs)) + " from " + str(len(args.input)) + " input files", file=sys.stderr)


    # Load reference and add labels
    ref = bed12.loadbed(args.reference, False, False)
    print ("# ref entries: " + str(len(ref)), file=sys.stderr)

    in_juncs = 0
    out_juncs = 0
    for i in range(0, len(bed)):
        b = bed[i]
        if b in ref:
            print("1")
        else:
            print("0")

main()