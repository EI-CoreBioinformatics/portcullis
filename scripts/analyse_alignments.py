#!/usr/bin/env python3

import os
from os.path import basename
import readline
import rpy2
import rpy2.robjects
from collections import OrderedDict, Counter
from rpy2.robjects.packages import importr
import itertools
import argparse
import bed12
import performance


def main():
    
    parser = argparse.ArgumentParser("Script to create the Venn Plots from BED files")
    parser.add_argument("input", nargs="+", help="The BED files to analyse")
    parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
    parser.add_argument("-o", "--output", required=True, help="The output prefix")
    args = parser.parse_args()

    ref_bed = bed12.loadbed(args.reference, False, False)
    print ("Loaded Reference BED file.  # junctions: " + str(len(ref_bed)))

    # Load all bed files
    bed_data = {}
    aligners = set()
    reads = set()
    #junc_analysers = set()
    for bed_path in args.input:
        bed_file = os.path.split(bed_path)[1]
        bed_base = os.path.splitext(bed_file)[0]
        bed_data[bed_base] = bed12.loadbed(bed_path, False, False)
        parts = bed_base.split('-')
        aligners.add(parts[0])
        reads.add(parts[1])
        #junc_analysers.add(parts[2])
        print ("Loaded: " + bed_file + "; # junctions: " + str(len(bed_data[bed_base])))

    print ("Found these aligners: " + ', '.join(aligners))
    print ("Found these reads: " + ', '.join(reads))
    #print ("Found these junction analysis tools: " + ', '.join(junc_analysers))


    # Build table
    tab = list()
    for a in aligners:
        for r in reads:
            p = performance.PEntry()
            p.aligner = a
            p.input = r
            p.tp = len(ref_bed & bed_data[a + "-" + r])
            p.fp = len(bed_data[a + "-" + r] - ref_bed)
            p.fn = len(ref_bed - bed_data[a + "-" + r])

            tab.append(r + "\t" + a + "\t" + p.__str__())

    # Output table to disk
    with open(args.output + "-align_reads.tab", "w") as tab_out:
        print("Dataset\tAligner\t" + performance.PEntry.header(), file=tab_out)
        for p in tab:
            print(p, file=tab_out)


    # Create Venns
    cols = rpy2.robjects.vectors.StrVector( ["lightblue", "purple", "green",
                                             "orange", "red"])

    r = rpy2.robjects.r  # Start the R thread
    base = importr("base")
    venn = importr("VennDiagram")
    grdevices = importr("grDevices")


    for r in reads:

        categories = list()
        categories.append("Reference")

        sets = list()
        sets.append(ref_bed)

        nums = dict()
        nums["area1"] = len(ref_bed)
        i=2
        for a in sorted(aligners):
            s = bed_data[a + "-" + r]
            sets.append(s)
            categories.append(a)
            nums["area{0}".format(i)] = len(s)
            i+=1

        for num_combs in range(2,6):
            for comb in itertools.combinations(range(1,6), num_combs):
                index = "".join([str(x) for x in comb])
                curr_sets = [sets[num-1] for num in comb]
                nums["n{0}".format(index)] = len(set.intersection(*curr_sets))

        grdevices.tiff(args.output + "-" + r + ".venn.tiff", width=960, height=960)
        venn.draw_quintuple_venn(height=5000,
                             width=5000,
                             # This will be in alphabetical order X(
                             fill=cols,
                             category=rpy2.robjects.vectors.StrVector(categories),
                             margin=0.2,
                             cat_dist=rpy2.robjects.vectors.FloatVector([0.25, 0.3, 0.25, 0.25, 0.25]),
                             cat_cex=3,
                             cat_col=rpy2.robjects.vectors.StrVector(["darkblue",
                                                                      "purple",
                                                                      "darkgreen",
                                                                      "darkorange",
                                                                      "darkred"]),
                             cex=2,
                             main="Comparison on junctions found by alignment tools",
                             main_col="black",
                             main_cex=8,
                             sub="" + r + " dataset",
                             sub_col="black",
                             sub_cex=5,
                             **nums)
    grdevices.dev_off()


main()
