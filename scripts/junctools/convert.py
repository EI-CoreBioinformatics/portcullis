#!/usr/bin/env python3

"""
This python script is intended to convert between various junction file formats.
"""

import argparse
import sys


from junctools.junction import *

__author__ = "Dan Mapleson"
__copyright__ = "Copyright 2016, Portcullis"
__credits__ = ["Dan Mapleson", "Luca Venturini", "David Swarbreck"]
__license__ = "GPLv3"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"



def decstart(junctions):
    for j in junctions:
        j.start -= 1


def tab2egff(args):
    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            l = line.strip()
            if not l == "":
                print(TabJunction().parse_line(l).toExonGFF(source=args.source))

def tab2igff(args):
    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            l = line.strip()
            if not l == "":
                print(TabJunction().parse_line(l).toIntronGFF(source=args.source))


def ebed2ibed(args):
    with open(args.input) as f:
        # Skip header
        f.readline()
        print("track name=\"junctions\"")

        for line in f:
            print(b = BedJunction(use_strand=True).parse_line(line).toIntronStyle())

def bed2ibed6(args):
    with open(args.input) as f:
        print(Bed6Junction.file_header())
        for line in f:
            b = BedJunction(use_strand=True, tophat=True).parse_line(line, fullparse=True)
            if not b == None:
                print(b.toIntronStyle(bed6=True))

def tbed2ebed(args):
    with open(args.input) as f:
        print(BedJunction.file_header())
        for line in f:
            b = BedJunction(use_strand=True, tophat=True).parse_line(line, fullparse=True)
            if not b == None:
                print(b)

def tbed2ibed(args):
    with open(args.input) as f:
        print(BedJunction.file_header())
        for line in f:
            b = BedJunction(use_strand=True, tophat=True).parse_line(line, fullparse=True)
            if not b == None:
                print(b.toIntronStyle())



def star2ibed(args):
    junctions = []

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            j = BedJunction()
            j.seq = words[0]
            j.start = int(words[1]) - 1
            j.end = int(words[2])
            j.left = j.start
            j.right = j.end
            j.strand = "+" if words[3] == "1" else "-" if words[3] == "2" else "."
            j.cov = int(words[6])
            junctions.append(j)

        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)


def fs2ibed(args):
    with open(args.input) as f:
        content = f.readline()

        print("track name=\"junctions\"")

        index = 0;
        for line in f:
            words = line.split()

            dif = int(words[2]) - int(words[1])

            print(words[0] + "\t" + str(words[1]) + "\t" + str(words[2]) + "\tjunc_" + str(index) + "\t0\t.\t" + str(
                words[1]) + "\t" + str(words[2]) + "\t255,0,0\t2\t1,1\t0," + str(dif))
            index += 1


def ts2ibed(args):
    junctions = []

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            j = BedJunction()
            j.seq = words[0]
            j.start = int(words[1]) - 1
            j.end = int(words[2]) - 1
            j.left = j.start
            j.right = j.end
            j.strand = "."
            j.cov = int(words[4])
            junctions.append(j)

        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)


def sp2ibed(args):
    junctions = []

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            parts1 = words[0].split(":")
            parts2 = parts1[1].split("_")

            j = BedJunction()
            j.seq = parts1[0]
            j.start = int(parts2[0]) - 1
            j.end = int(parts2[1])
            j.left = j.start
            j.right = j.end
            j.strand = parts1[2]
            j.cov = int(words[9])
            junctions.append(j)

        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)


def ss2ibed(args):
    junctions = list()

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            j = BedJunction()
            j.seq = words[0]
            j.start = int(words[1])
            j.end = int(words[2]) - 1
            j.left = j.start
            j.right = j.end
            j.strand = "."
            j.cov = int(words[4])
            junctions.append(j)

        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)

def ms2ibed(args):
    junctions = list()

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            j = Junction()
            j.seq = words[0]
            j.start = int(words[1])
            j.end = int(words[2]) - 1
            j.left = j.start
            j.right = j.end
            j.strand = words[5]
            j.cov = int(words[4])
            junctions.append(j)

        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)

def ebed2hisat(args):
    with open(args.input) as f:
        for line in f:
            b = BedJunction().parse_line(line)
            if not b == None:
                print("\t".join([b.refseq, str(b.start-1), str(b.end+1), b.strand]))

def gtf2ibed(args):
    with open(args.input) as f:
        junctions = [];
        junction_set = {}

        curr_transcript_seq = ""
        curr_transcript_start = 0
        curr_transcript_end = 0
        curr_transcript_strand = ""

        last_exon_seq = ""
        last_exon_start = 0
        last_exon_end = 0
        last_exon_strand = ""

        index = 0;
        for line in f:

            if not line.startswith("#"):
                words = line.split()

                start = int(words[3])
                end = int(words[4])

                dif = end - start

                if words[2] == "transcript":

                    curr_transcript_seq = words[0]
                    curr_transcript_start = start
                    curr_transcript_end = end
                    curr_transcript_strand = words[6]

                    # Wipe last exon
                    last_exon_seq = ""
                    last_exon_start = 0
                    last_exon_end = 0
                    last_exon_strand = ""

                elif words[2] == "exon":

                    # Double check we are still in the current transcript and there has already been an exon in this transcript
                    if curr_transcript_end >= last_exon_end and curr_transcript_seq == last_exon_seq and curr_transcript_strand == last_exon_strand:

                        j = BedJunction()
                        j.seq = words[0]
                        j.start = last_exon_end + 1
                        j.end = start - 1
                        j.strand = words[6]

                        if j.key() not in junction_set:
                            junctions.append(j)
                            junction_set[j.key()] = j

                    last_exon_seq = words[0]
                    last_exon_start = start
                    last_exon_end = end
                    last_exon_strand = words[6]

        decstart(junctions)
        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)


def gff2ibed(args):
    with open(sys.argv[1]) as f:
        content = f.readline()

        junctions_set = {}
        junctions = []

        index = 0;
        for line in f:

            if not line.startswith("#"):
                words = line.split()

                if words[2] == "intron":

                    start = int(words[3]) - 1
                    end = int(words[4])

                    dif = end - start
                    key = words[0] + "_" + str(start) + "_" + str(end) + "_" + words[6]
                    line = words[0] + "\t" + str(start) + "\t" + str(end) + "\tjunc_" + str(index) + "\t0\t" + words[
                        6] + "\t" + str(start) + "\t" + str(end) + "\t255,0,0\t2\t0,0\t0," + str(dif)

                    j = Junction()
                    j.seq = words[0]
                    j.start = start
                    j.end = end
                    j.strand = words[6]

                    if key not in junctions_set:
                        junctions_set[key] = line
                        junctions.append(j)
                        index += 1

        Junction.sort(junctions)
        Junction.reindex(junctions)

        print(BedJunction.file_header())
        for j in junctions:
            print(j)





