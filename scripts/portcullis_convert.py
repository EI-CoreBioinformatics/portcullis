#!/usr/bin/env python3

"""
This python script is intended to convert between various junction file formats.
"""

import sys
import argparse

import bed12


__author__ = "Dan Mapleson"
__copyright__ = "Copyright 2016, Portcullis"
__credits__ = ["Dan Mapleson", "Luca Venturini", "David Swarbreck"]
__license__ = "GPLv3"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"


class Junction:
    def __init__(self):
        self.data = []
        self.seq = ""
        self.start = 0
        self.end = 0
        self.strand = "."
        self.cov = 0

    def display(self):
        print(self.seq + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + self.strand + "\t" + str(self.cov))

    def key(self):
        return (self.seq, self.start, self.end, self.strand)


def decstart(junctions):
    for j in junctions:
        j.start -= 1


def sortandprint(junctions):
    # Sort by
    junctions.sort(key=lambda x: x.strand)
    junctions.sort(key=lambda x: x.end)
    junctions.sort(key=lambda x: x.start)
    junctions.sort(key=lambda x: x.seq)

    print("track name=\"junctions\"")

    index = 0;
    for j in junctions:
        print(j.seq + "\t" + str(j.start) + "\t" + str(j.end) + "\tjunc_" + str(index) + "\t" + str(
            j.cov) + "\t" + j.strand + "\t" + str(j.start) + "\t" + str(j.end) + "\t255,0,0\t2\t0,0\t0,0")

        index += 1


def ebed2ibed(args):
    with open(args.input) as f:
        # Skip header
        f.readline()
        print("track name=\"junctions\"")

        for line in f:
            b = bed12.BedEntry.create_from_line(line, use_strand=True)
            print(b.toIntronStyle())

def ebed2ibed6(args):
    with open(args.input) as f:
        # Skip header
        f.readline()
        print("track name=\"junctions\"")

        for line in f:
            b = bed12.BedEntry.create_from_line(line, use_strand=True)
            print(b.toIntronStyle(bed6=True))

def ibed2ibed6(args):
    with open(args.input) as f:
        # Skip header
        f.readline()
        print("track name=\"junctions\"")

        for line in f:

            parts = line.split("\t")
            print("\t".join(parts[0:6]))


def tbed2ebed(args):
    with open(args.input) as f:
        # Skip header
        f.readline()
        print("track name=\"junctions\"")

        for line in f:
            print(bed12.BedEntry.create_from_line(line, use_strand=True, tophat=True))

def tbed2ibed(args):
    with open(args.input) as f:
        # Skip header
        f.readline()
        print("track name=\"junctions\"")

        for line in f:
            b = bed12.BedEntry.create_from_line(line, use_strand=True, tophat=True)
            print(b.toIntronStyle())



def star2ibed(args):
    junctions = []

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            j = Junction()
            j.seq = words[0]
            j.start = int(words[1]) - 1
            j.end = int(words[2])
            j.strand = "+" if words[3] == "1" else "-" if words[3] == "2" else "."
            j.cov = int(words[6])
            junctions.append(j)

        sortandprint(junctions)


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

            j = Junction()
            j.seq = words[0]
            j.start = int(words[1]) - 1
            j.end = int(words[2]) - 1
            j.strand = "."
            j.cov = int(words[4])
            junctions.append(j)

        sortandprint(junctions)


def sp2ibed(args):
    junctions = []

    with open(args.input) as f:
        # Skip header
        f.readline()

        for line in f:
            words = line.split("\t")

            parts1 = words[0].split(":")
            parts2 = parts1[1].split("_")

            j = Junction()
            j.seq = parts1[0]
            j.start = int(parts2[0]) - 1
            j.end = int(parts2[1])
            j.strand = parts1[2]
            j.cov = int(words[9])
            junctions.append(j)

        sortandprint(junctions)


def ss2ibed(args):
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
            j.strand = "."
            j.cov = int(words[4])
            junctions.append(j)

        sortandprint(junctions)

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
            j.strand = words[5]
            j.cov = int(words[4])
            junctions.append(j)

        sortandprint(junctions)

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

                        j = Junction()
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
        sortandprint(junctions)


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

        sortandprint(junctions)


def main():
    call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        """This script contains a set of options for converting between various
    splice junction format files.""")
    subparsers = parser.add_subparsers(
        title="Conversion options")

    ebed2ibed_parser = subparsers.add_parser("ebed2ibed",
                                             help="Converts a portcullis BED file (containing exon anchors - as produced by the main executable) to a pure intron-based BED 12 file (no exon anchors).")

    ebed2ibed_parser.add_argument("input", help="The portcullis exon-based BED file with anchors to convert")
    ebed2ibed_parser.set_defaults(func=ebed2ibed)

    ebed2ibed6_parser = subparsers.add_parser("ebed2ibed6",
                                             help="Converts a portcullis BED file (containing exon anchors - as produced by the main executable) to a pure intron-based BED 6 file (no exon anchors).")
    ebed2ibed6_parser.add_argument("input", help="The portcullis exon-based BED file with anchors to convert")
    ebed2ibed6_parser.set_defaults(func=ebed2ibed6)

    ibed2ibed6_parser = subparsers.add_parser("ibed2ibed6",
                                             help="Converts a portcullis intron-based BED 12 file (12 column format) to BED 6 file.  This tools simply removes the last 6 columns from the BED file.")
    ibed2ibed6_parser.add_argument("input", help="The intron-based (no anchors) BED 12 file")
    ibed2ibed6_parser.set_defaults(func=ibed2ibed6)


    ebed2ibed_parser.add_argument("input", help="The portcullis BED file to convert")
    ebed2ibed_parser.set_defaults(func=ebed2ibed)

    tbed2ebed_parser = subparsers.add_parser("tbed2ebed",
                                             help="Converts a tophat BED file (containing exon anchors - junction sites defined by block start+size) to a portcullis exon-based BED file (exon anchors - thickstart+end define junction sites).")

    tbed2ebed_parser.add_argument("input", help="The tophat BED file to convert")
    tbed2ebed_parser.set_defaults(func=tbed2ebed)

    tbed2ibed_parser = subparsers.add_parser("tbed2ibed",
                                             help="Converts a tophat BED file (containing exon anchors - junction sites defined by block start+size) to a portcullis intron-based BED file (no anchors).")

    tbed2ibed_parser.add_argument("input", help="The tophat BED file to convert")
    tbed2ibed_parser.set_defaults(func=tbed2ibed)

    star2ibed_parser = subparsers.add_parser("star2ibed",
                                             help="Converts a STAR junction file to a portcullis (intron-based - no anchors) BED file.")

    star2ibed_parser.add_argument("input", help="The STAR tab delimited junction file to convert")
    star2ibed_parser.set_defaults(func=star2ibed)

    fs2ibed_parser = subparsers.add_parser("finesplice2ibed",
                                           help="Converts a finesplice junction file to a portcullis (intron-based - no anchors) BED file.")

    fs2ibed_parser.add_argument("input", help="The finesplice tab delimited junction file to convert")
    fs2ibed_parser.set_defaults(func=fs2ibed)

    ts2ibed_parser = subparsers.add_parser("truesight2ibed",
                                           help="Converts a truesight junction file to a portcullis (intron-based - no anchors) BED file.")

    ts2ibed_parser.add_argument("input", help="The truesight tab delimited junction file to convert")
    ts2ibed_parser.set_defaults(func=ts2ibed)

    sp2ibed_parser = subparsers.add_parser("spanki2ibed",
                                           help="Converts a SPANKI junction file to a portcullis (intron-based - no anchors) BED file.")

    sp2ibed_parser.add_argument("input", help="The SPANKI tab delimited junction file to convert")
    sp2ibed_parser.set_defaults(func=sp2ibed)

    ss2ibed_parser = subparsers.add_parser("soapsplice2ibed",
                                           help="Converts a soapsplice junction file to a portcullis (intron-based - no anchors) BED file.")

    ss2ibed_parser.add_argument("input", help="The soapsplice tab delimited junction file to convert")
    ss2ibed_parser.set_defaults(func=ss2ibed)
    
    ms2ibed_parser = subparsers.add_parser("mapsplice2ibed",
                                           help="Converts a mapsplice junction file to a portcullis (intron-based - no anchors) BED file.")

    ms2ibed_parser.add_argument("input", help="The mpsplice tab delimited junction file to convert")
    ms2ibed_parser.set_defaults(func=ms2ibed)

    gtf2ibed_parser = subparsers.add_parser("gtf2ibed",
                                            help="Converts a GTF file containing genes and transcripts in to a portcullis intron-based BED file containing junctions derived from the GTF.")

    gtf2ibed_parser.add_argument("input", help="The GTF file from which to extract junctions")
    gtf2ibed_parser.set_defaults(func=gtf2ibed)

    gff2ibed_parser = subparsers.add_parser("gff2ibed",
                                            help="Converts a GFF file containing introns to a portcullis intron-based BED file. Check that GFF file has been marked up with introns first!  If not, you can use genome tools for this.")

    gff2ibed_parser.add_argument("input", help="The GFF file containing introns to convert into BED")
    gff2ibed_parser.set_defaults(func=gff2ibed)

    args = parser.parse_args(call_args)
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
