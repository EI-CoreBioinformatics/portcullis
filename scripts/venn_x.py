import rpy2.robjects
# import sys
import csv
import os
import rpy2.robjects
from collections import OrderedDict, Counter
from rpy2.robjects.packages import importr
import itertools
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

def main():
    
    parser = argparse.ArgumentParser("Script to create the Venn Plots from BED files")
    parser.add_argument("input", nargs='+', help="The input BED files")
    args = parser.parse_args()

    bed_list = ()
    for b in args.input:
        bed_list.append(loadbed(b, false, false))

    print(str(len(bed_list)))

    '''
    r = rpy2.robjects.r  # Start the R thread                                                                                       
    base = importr("base")
    venn = importr("VennDiagram")
    grdevices = importr("grDevices")
    # corrs = {1: "class", 2: "cufflinks", 3: "stringtie", 4: "trinity", 5: "mikado"}
    corrs = dict((x+1, list(sets.keys())[x]) for x in range(len(sets.keys())))
    nums = dict()

    for num in corrs:
        cat = corrs[num]
        nums["area{0}".format(num)] = len(sets[cat])
        total.update(list(sets[cat]))
        if num < 5:
            total_wo_mikado.update(list(sets[cat]))
        print(cat.capitalize(), nums["area{0}".format(num)])

    print("Total", len(set.union(*sets.values())))
    print("Total w/o Mikado", len(set.union(*[sets[x] for x in sets if x != "mikado"])))
    #
    counts = list(total.values())
    print("## With Mikado")
    for num in reversed(range(6)):
        tot = counts.count(num)
        cum_tot = sum( counts.count(x) for x in range(num+1, 6)) + tot
        print("Genes {0} by {1} methods ({2} cumulative)".format(args.type, num,
                                                                 cum_tot), tot)
    print("")
    print("## Without Mikado")
    counts = list(total_wo_mikado.values())
    for num in reversed(range(5)):
        tot = counts.count(num)
        cum_tot = sum( counts.count(x) for x in range(num+1, 5)) + tot
        print("Genes {0} by {1} methods ({2} cumulative)".format(args.type, num,
                                                                 cum_tot), tot)

    for num_combs in range(2,6):
        for comb in itertools.combinations(range(1,6), num_combs):
            index = "".join([str(x) for x in comb])
            curr_sets = [sets[corrs[num]] for num in comb]
            nums["n{0}".format(index)] = len(set.intersection(*curr_sets))

    cols = rpy2.robjects.vectors.StrVector( ["lightblue", "purple", "green",
                                             "orange", "red"])

    grdevices.tiff(os.path.join(args.out_folder, "quintuple_{0}.new.tiff".format(args.type)),
                   width=960, height=960)

    venn.draw_quintuple_venn(height=5000,
                             width=5000,
                             # This will be in alphabetical order X(
                             fill=cols,
                             category=rpy2.robjects.vectors.StrVector([x.capitalize() for x in sets.keys()]),
                             margin=0.2,
                             cat_dist=rpy2.robjects.vectors.FloatVector([0.25, 0.3, 0.25, 0.25, 0.25]),
                             cat_cex=3,
                             cat_col=rpy2.robjects.vectors.StrVector(["darkblue",
                                                                      "purple",
                                                                      "darkgreen",
                                                                      "darkorange",
                                                                      "darkred"]),
                             cex=2,
                             **nums)
    grdevices.dev_off()

    grdevices.tiff(os.path.join(args.out_folder,
                                "quadruple_{0}.new.tiff".format(args.type)),
                   width=960, height=960)
    cols = rpy2.robjects.vectors.StrVector(["lightblue", "purple", "green",
                                            "chocolate2"])

    nums = dict((x, nums[x]) for x in nums.keys() if "5" not in x)

    title_vector = []
    denominator = len(set.union(*[sets[x] for x in sets if x != "mikado"]))
    print(nums.keys())
    title_vector.append("Class\n{0:,} genes ({1}%)".format(nums["area1"],
                                                         round(100*nums["area1"]/denominator,1)))
    title_vector.append("Cufflinks\n{0:,} genes ({1}%)".format(nums["area2"],
                                                             round(100*nums["area2"]/denominator,1)))
    title_vector.append("Stringtie\n{0:,} genes ({1}%)".format(nums["area3"],
                                                             round(100*nums["area3"]/denominator,1)))
    title_vector.append("Trinity\n{0:,} genes ({1}%)".format(nums["area4"],
                                                           round(100*nums["area4"]/denominator,1)))
    
    venn.draw_quad_venn(height=5000,
                        width=5000,
                        # This will be in alphabetical order X(
                        fill=cols,
                        # category = rpy2.robjects.vectors.StrVector(["Class", "Cufflinks",
                        # "Stringtie", "Trinity"]),
                        category=rpy2.robjects.vectors.StrVector(title_vector),
                        margin=0.15,
                        cat_dist=rpy2.robjects.vectors.FloatVector([0.32, 0.3, 0.2, 0.25]),
                        cat_cex=2.8,
                        cat_col=rpy2.robjects.vectors.StrVector(
                            ["darkblue", "purple", "darkgreen", "chocolate4"]),
                        cex = 2.2,
                        **nums)
    grdevices.dev_off()
    '''
main()
