#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

import sys
import argparse
import copy
import collections
# from Mikado.parsers.bed12 import Bed12Parser, BED12
import logging
# from Mikado.utilities.log_utils import create_default_logger
from collections import namedtuple
from bed12 import *

def calcOp(op, newval, currval) :
    if op == "max":
        return max(newval, currval)
    elif op == "total":
        return newval + currval
    else:
        # Default to total
        return newval + currval

# @profile
def main():
    parser=argparse.ArgumentParser("Script to merge multiple BED files, using target sequence, and thick_start and thick_end values to generate set keys.  (Only entering a single file here essentially copies the file).  Merged BED entries contain the maximum score found across all input files.")
    parser.add_argument("input", nargs="+", help="List of BED files to intersect")
    parser.add_argument("-m", "--min_entry", type=int, default="1", help="Minimum number of files the entry is require to be in.  0 means entry must be present in all files, i.e. true intersection.  1 means a union of all input files")
    parser.add_argument("--operator", default="total", help="Operator to use for calculating the score in the merged BAM file.  Options: [max, total]")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument("-l", "--log", default=None, type=argparse.FileType("w"))
    parser.add_argument("-s", "--ignore_strand", action='store_true', default=False, help="Ignore strand information for calculating match in bed files")
    parser.add_argument("-p", "--prefix", default="Junc_intersected", help="Prefix to apply to name column in BED output file")
    args=parser.parse_args()

    min_entry = args.min_entry if args.min_entry > 0 else len(args.input)
    logger = logging.getLogger("merger")
    formatter = logging.Formatter(
            "{asctime} - {name} - {filename}:{lineno} - {levelname} - {funcName} - {message}",
            style="{"
            )
    logger.setLevel(logging.INFO)
    if args.log is not None:
        args.log.close()
        args.log = args.log.name
        handler = logging.FileHandler(args.log)
    else:
        handler = logging.StreamHandler()
        
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    logger.info("Requested %d input BED files to be intersected", len(args.input))
    logger.info("Each entry must be found in at least %d files to be in output", min_entry)

    bed_merge = collections.defaultdict(list)
    
    for bf in args.input:
        counter = 0
        found = set()
        with open(bf) as bed_file:
            counter = 0
            found = set()
            for line in bed_file:
                entry = BedEntry.create_from_line(line, not args.ignore_strand, False)
                if entry is None:
                    # Skipping header
                    continue
                counter += 1
                found.add(entry.key)
                bed_merge[entry.key].append(line)
                continue
            pass
                
        if len(found) < counter:
            logger.warning("File %s had duplicated entries (%d lines, %d unique entries)",
                                       bf, counter, len(found))
        logger.info("Input file %s contains %d entries.", bf, len(bed_merge))

    del found
    logger.info("Found a total of %d distinct entries.", len(bed_merge))

    # bed_out = list()
    i = 1
    with open(args.output, "wt") as out:
        description = "Merge of multiple BED files.  Min_Entry: {0}. Score_op: {1}".format(min_entry, args.operator.upper())
        print("track name=\"junctions\" description=\"{0}\"".format(description), file=out)
        for b in sorted(bed_merge):
            if len(bed_merge[b]) >= min_entry:
                bo = BedEntry.create_from_line(bed_merge[b][0], not args.ignore_strand, False)
                bo.name = "{prefix}_{i}".format(prefix=args.prefix, i=i)
                i += 1
                for r in bed_merge[b][1:]:
                    bo.score = calcOp(args.operator,
                                      BedEntry.create_from_line(r, not args.ignore_strand, False).score,
                                      bo.score)
                # for num in (8,10,11):
                #     line[num] = ",".join([str(_) for _ in line[num]])
                print(bo, file=out)

    logger.info("Filtered out %d entries", len(bed_merge) - i + 1)
    logger.info("Output file %s contains %d entries.", args.output, counter - 1)
    if args.log is not None:
        handler.close()
        
if __name__ == '__main__':
    main()
