#!/usr/bin/env python3

__author__ = 'maplesod'

from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
import sys
import argparse

def loadbed(filepath, usestrand, tophat) :
    with open(filepath) as f:
        index = 0
        items = set()
        for line in f:
            if index > 0:
                words = line.split()
                overhang = words[10]
                overhang_parts = overhang.split(",")
                lo = int(overhang_parts[0])
                ro = int(overhang_parts[1])
                chr = words[0]
                start = str(int(words[6]) + lo) if tophat else words[6]
                end = str(int(words[7]) - ro) if tophat else words[7]
                strand = words[5]
                key = chr + "_" + start + "_" + end
                if usestrand:
                    key += "_" + strand
                items.add(key)
            index += 1
    if len(items) != index - 1 :
        print ("non unique items in bed file " + filepath)
    return items

def filtbed(infile, outfile, keeplist, usestrand) :
    o = open(outfile, 'w')
    with open(infile) as f:
        index=0
        for line in f:
            if index > 0:
                words = line.split()
                chr = words[0]
                start = words[6]
                end = words[7]
                strand = words[5]
                key = chr + "_" + start + "_" + end
                if usestrand:
                    key += "_" + strand
                if key in keeplist:
                    print >> o, line,
            index += 1
    o.close()
    return

def filttab(tabfile, bedfile, filttabfile) :

    items = set()
    with open(bedfile) as b:
        index = 0
        for line in b:
            if index > 0:
                words = line.split()
                chr = words[0]
                start = words[6]
                end = words[7]
                strand = words[5]
                key = chr + "_" + start + "_" + end + "_" + strand
                items.add(key)
            index +=1

    print (len(items))

    o = open(filttabfile, 'w')
    with open(tabfile) as f:
        index=0
        for line in f:
            if index != 0 and not line.isspace():
                words = line.split()
                chr = words[2]
                start = words[4]
                end = words[5]
                strand = words[11]
                key = chr + "_" + start + "_" + end + "_" + strand
                if key in items:
                    print >> o, line,
            else:
                print >> o, line,
            index += 1
    o.close()
    return

    
parser=argparse.ArgumentParser("Script to analyse the results of portcullis when run on data from a model organism (i.e. it has a high quality annotated genome.")
parser.add_argument("-r", "--reference", required=True,
                    help="Reference bed file.")
parser.add_argument("-p", "--potential", required=True,
                    help="All potential junctions in bed format.  Take output from `portcullis junc`")#
parser.add_argument("-f", "--filtered", required=True,
                    help="Junctions after filtering in bed format.")
parser.add_argument("-o", "--out", required=True,
                    help="Output for Venn diagrams")
parser.add_argument("-t", "--tophat", action='store_true', default=False,
                    help="Filtered junctions is from tophat, compensate for maximal overhangs")
parser.add_argument("-s", "--ignore_strand", action='store_true', default=False,
                    help="Use strand information in bed files")
args=parser.parse_args()

print ("Reference: " + args.reference)
print ("Potential aligned junctions: " + args.potential)
print ("Filtered junctions: " + args.filtered)
print ("Output PNG: " + args.out)

trueset = loadbed(args.reference, not args.ignore_strand, False)
alignedset = loadbed(args.potential, not args.ignore_strand, False)
filteredset = loadbed(args.filtered, not args.ignore_strand, args.tophat)

plt.figure(figsize=(10,10))
#v = venn3([trueset, alignedset, filteredset], set_labels = (sys.argv[1], sys.argv[2], sys.argv[3]))
v = venn3_unweighted([trueset, alignedset, filteredset], set_labels = (args.reference, args.potential, args.filtered))
#c = venn3_circles([trueset, alignedset, filteredset], linestyle='dashed')
#c[0].set_lw(1.0)
#c[0].set_ls('dotted')

# Work out the scores
true_juncs = trueset.intersection(alignedset)
false_juncs = alignedset.difference(trueset)
true_positives = filteredset.intersection(true_juncs)
false_positives = filteredset.intersection(false_juncs)
unrecoverable_juncs = trueset.difference(alignedset)

true_negatives = false_juncs.difference(filteredset)
false_negatives = true_juncs.difference(filteredset)

tp = len(true_positives)
fp = len(false_positives)
tn = len(true_negatives)
fn = len(false_negatives)

print
print ("Junctions in reference: " + str(len(trueset)))
print ("Potential true junctions: " + str(len(true_juncs)))
print ("Unrecoverable junctions (not in aligned data): " + str(len(unrecoverable_juncs)))
print ("False junctions (not in reference): " + str(len(false_juncs)))
print

print ("True positives: " + str(tp))
print ("False positives: " + str(fp))
print ("True negatives: " + str(tn))
print ("False negatives: " + str(fn))
print

sen = float(tp) / float(tp + fn)
spc = float(tn) / float(tn + fp)
prc = float(tp) / float(tp + fp)

print ("Sensitivity: " + str(sen))
print ("Specificity: " + str(spc))
print ("Precision:" + str(prc))
print


acc = float(tp + tn) / float(tp + fp + tn + fn)
f1 = float(2 * tp) / float(2 * tp + fp + fn)

print ("Accuracy: " + str(acc))
print ("F1 Score: " + str(f1))

plt.title(args.out)
plt.savefig(args.out)

sys.exit(0)