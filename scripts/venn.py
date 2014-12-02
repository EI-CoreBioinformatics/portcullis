__author__ = 'maplesod'

from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
import sys
def loadbed(filepath, usestrand) :
    with open(filepath) as f:
        index = 0
        items = set()
        for line in f:
            if index != 0:
                words = line.split()
                chr = words[0]
                start = words[6]
                end = words[7]
                strand = words[5]
                key = chr + "_" + start + "_" + end
                if usestrand:
                    key += "_" + strand
                items.add(key)
            index += 1
    if len(items) != index :
        print "non unique items in bed file " + filepath
    return items

def filtbed(infile, outfile, keeplist, usestrand) :
    o = open(outfile, 'w')
    with open(infile) as f:
        index=0
        for line in f:
            if index != 0:
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
            if index != 0:
                words = line.split()
                chr = words[0]
                start = words[6]
                end = words[7]
                strand = words[5]
                key = chr + "_" + start + "_" + end + "_" + strand
                items.add(key)
            index +=1

    print len(items)

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



trueset = loadbed(sys.argv[1], 0)
alignedset = loadbed(sys.argv[2], 0)
filteredset = loadbed(sys.argv[3], 0)

plt.figure(figsize=(10,10))
#v = venn3([trueset, alignedset, filteredset], set_labels = (sys.argv[1], sys.argv[2], sys.argv[3]))
v = venn3_unweighted([trueset, alignedset, filteredset], set_labels = (sys.argv[1], sys.argv[2], sys.argv[3]))
#c = venn3_circles([trueset, alignedset, filteredset], linestyle='dashed')
#c[0].set_lw(1.0)
#c[0].set_ls('dotted')
plt.title(sys.argv[4])
plt.savefig(sys.argv[4])

