#!/usr/bin/env python3

import itertools
import argparse
import bed12
from performance import Performance
from matplotlib_venn import venn2
from pylab import figure


def main():
    parser = argparse.ArgumentParser("Script to quantify potential alternative splicing events within a single sample")
    parser.add_argument("input", nargs="+", help="The BED file(s) to analyse")
    args = parser.parse_args()

    print("File\tjuncs\tdupstart\tdupend\ttotaldup\tmeandupcov")

    for bf in args.input:
        bed_data = bed12.loadbed(bf, False, False)

        start = set()
        end = set()

        dup_start = 0
        dup_end = 0
        dup_sum = 0
        dup_count = 0

        for b in bed_data:
            if b.thick_start in start:
                dup_start += 1
                dup_sum += b.score
                dup_count += 1
            else:
                start.add(b.thick_start)
            if b.thick_end in end:
                dup_end += 1
                dup_sum += b.score
                dup_count += 1
            else:
                end.add(b.thick_end)


        print(bf, "\t", len(bed_data), "\t", dup_start, "\t", dup_end, "\t", dup_count, "\t", float(dup_sum) / float(dup_count))


main()
