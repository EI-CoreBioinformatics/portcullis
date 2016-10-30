#!/usr/bin/env python3

import argparse
import sys

from junctools.junction import JuncFactory, Junction
from junctools import merge


def convert(args):
    in_type = JuncFactory[args.input_format.upper()]
    out_type = JuncFactory[args.output_format.upper()]

    if out_type == JuncFactory.TAB:
        raise "Can't output to TAB format"

    junction_set = {}
    junctions = []
    with open(args.input) as f:
        for line in f:
            j = JuncFactory.create_from_enum(in_type, use_strand=not args.ignore_strand).parse_line(line)

            if j:

                if args.dedup:
                    if j.key not in junction_set:
                        junctions.append(j)
                else:
                    junctions.append(j)

    if args.sort:
        Junction.sort(junctions)

    if args.reindex:
        Junction.reindex(junctions, prefix=args.prefix, start=args.index_start)

    o = sys.stdout
    if args.output != sys.stdout:
        o = open(args.output, mode='w')

    header = JuncFactory.create_from_enum(out_type).file_header()

    if header and header != "":
        print(header, file=o)

    for j in junctions:

        # Convert
        c = JuncFactory.create_from_enum(out_type, use_strand=not args.ignore_strand, junc_to_copy=j)

        if out_type.isBed():
            print(c, style=out_type, file=o)
        else:
            print(c, file=o)


def add_merge_options(parser):
    # parser = argparse.ArgumentParser("This script can merge output from several portcullis runs into one.\
    #		This supports BED12 and TAB formats that contain exon anchors.\
    #		Please note for TAB format that while anchors will be extended the metrics used to a merged junction will reflect\
    #		only the junction from the first sample and are therefore not recalculated based on the merged information")

    parser.add_argument("-m", "--min_entry", type=int, default="1",
                        help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.  1 means a union of all input files")
    parser.add_argument("--operator", default="total",
                        help="Operator to use for calculating the score in the merged file.  Options: [min, max, sum, mean]")
    parser.add_argument("-o", "--output", required=True, help="Merged output file")
    parser.add_argument("-p", "--prefix", default="junc_merged",
                        help="Prefix to apply to name column in BED output file")
    parser.add_argument("input", nargs="+", help="List of BED or TAB files to merge (must all be the same type)")

    parser.set_defaults(func=merge.merge)


def main():
    call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        """This script contains a number of tools for manipulating junction files.""",
        formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(
        title="Junction tools")

    convert_parser = subparsers.add_parser("convert",
                                           help="Converts junction files between various formats.")

    convert_parser.formatter_class = argparse.RawTextHelpFormatter

    convert_parser.add_argument("-if", "--input_format", required=True, help='''The format of the input file to convert.  Available options:
tab    = Portcullis tab output
bed    = BED format (We automatically determine if this is BED 6 or 12 format,
         and if it is intron, exon or tophat style)
star   = STAR style tab delimited format
hisat  = HISAT style tab delimited format
fs     = Finesplice style tab delimited format
ms     = Mapsplice style tab delimited format
ss     = Soapsplice style tab delimited format
ts     = Truesight style tab delimited format
spanki = SPANKI style tab delimited format
gtf    = Transcript assembly or gene model.  Note: output will only contain
         junctions
gff    = Transcript assembly or gene model containing introns to extract.
         Note: input must contain \"intron\" entries, and output will only
         contain these introns represented as junctions'''
                                )

    convert_parser.add_argument("-of", "--output_format", required=True, help='''The output format.  Available options:
ebed   = Exon-based BED12 format (Thick-start and end represent splice sites)
tbed   = Tophat style exon-based BED12 format (splice sites derived from blocks)
ibed   = Intron-based BED12 format
bed6   = BED6 format (BED6 files are intron-based)
hisat  = HISAT style tab delimited format
egff   = Exon-based junctions in GFF3 format, uses partial matches to indicate exon anchors
igff   = Intron-based junctions in GFF3 format'''
                                )

    convert_parser.add_argument("-o", "--output", default=sys.stdout, help="Output to this file.  By default we print to stdout.")

    convert_parser.add_argument("-is", "--ignore_strand", action='store_true', default=False,
                        help="Whether or not to ignore strand when creating a key for the junction")

    convert_parser.add_argument("-d", "--dedup", action='store_true', default=False,
                        help="Whether or not to remove duplicate junctions")

    convert_parser.add_argument("-s", "--sort", action='store_true', default=False,
                        help="Whether or not to sort the junctions")

    convert_parser.add_argument("-r", "--reindex", action='store_true', default=False,
                        help="Whether or not to reindex the output.  The index is applied after prefix.")

    convert_parser.add_argument("--index_start", type=int, default=0,
                        help="The starting index to apply if the user requested reindexing")

    convert_parser.add_argument("--prefix", default="junc_",
                        help="The prefix to apply to junction ids if the user requested reindexing")

    convert_parser.add_argument("--source", default="portcullis",
                        help="Only relevant if output is GFF format, use this option to set the source column in the GFF")

    convert_parser.add_argument("input", help="The input file to convert")

    convert_parser.set_defaults(func=convert)

    merge_parser = subparsers.add_parser("merge",
                                         help="Merges several junction files into a single junction file.")
    add_merge_options(merge_parser)

    args = parser.parse_args(call_args)
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
