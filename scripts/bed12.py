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



def filterbed(filepath, refset, mode, usestrand, tophat, outfile) :

    o = open(outfile, 'w')

    o.write("track name=\"junctions\"\n")

    index = 0
    with open(filepath) as i:

        i.readline()    # Skip header
        for line in i:

            key = makekey(line, usestrand, tophat)

            if mode == 0:
                if key in refset:
                    o.write(line)
            elif mode == 1:
                if key not in refset:
                    o.write(line)

            index += 1
    o.close()

def filtertab(filepath, outfile, bed, mode, usestrand) :

    o = open(outfile, 'w')

    index = 0;
    with open(filepath) as f:

        o.write(f.readline())

        for line in f:

            line = line.strip()
            if line != "":
                key = makekeyfromtab(line, usestrand)
                if mode == 0:
                    if key in bed:
                        o.write(line + "\n")
                elif mode == 1:
                    if key not in bed:
                        o.write(line + "\n")

    o.close()

