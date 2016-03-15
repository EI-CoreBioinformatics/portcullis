#!/usr/bin/env python3

import itertools
import argparse
import bed12
from matplotlib_venn import venn2
from pylab import figure

class PEntry:

    def __init__(self):
        self.tp = 0
        self.fp = 0
        self.fn = 0

    def __str__(self):
        return str(self.tp) + "\t" + str(self.fp) + "\t" + str(self.fn) + "\t" + format(self.calcRecall(), '.2f') \
                + "\t" + format(self.calcPrecision(), '.2f') + "\t" + format(self.calcF1(), '.2f')

    def calcPrecision(self):
        return 100.0 * float(self.tp) / float(self.tp + self.fp)

    def calcRecall(self):
        return 100.0 * float(self.tp) / float(self.tp + self.fn)

    def calcF1(self):
        prc = self.calcPrecision()
        rec = self.calcRecall()
        return 2.0 * (prc * rec) / (prc + rec)

    @staticmethod
    def header():
        return "TP\tFP\tFN\tRecall\tPrecision\tF1"

def main():
    
    parser = argparse.ArgumentParser("Script to compare bed file against reference bed")
    parser.add_argument("input", help="The BED file to analyse")
    parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
    parser.add_argument("-o", "--output", required=True, help="The output prefix")
    args = parser.parse_args()

    ref_bed = bed12.loadbed(args.reference, False, False)
    print ("Loaded Reference BED file.  # junctions: ", len(ref_bed))

    # Load all bed files
    bed_data = bed12.loadbed(args.input, False, False)
    print ("Loaded: " + args.input + "; # junctions: ", len(bed_data))

    # Build table
    tab = list()
    p = PEntry()
    p.tp = len(ref_bed & bed_data)
    p.fp = len(bed_data - ref_bed)
    p.fn = len(ref_bed - bed_data)

    print("Results:")
    print(PEntry.header())
    print(p)

    # Create Venns
    plt = figure(1, figsize=(6, 6))
    venn2(subsets=(p.fn, p.fp, p.tp), set_labels=(args.reference, args.input))
    plt.show()
    plt.savefig(args.output)

main()
