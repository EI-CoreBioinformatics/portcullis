#!/usr/bin/env python3

import argparse

from matplotlib_venn import venn2
from pylab import figure

from junctools.performance import Performance
from junctools.junction import *


def main():
	parser = argparse.ArgumentParser("Script to compare bed file against reference bed")
	parser.add_argument("input", nargs="+", help="The BED file to analyse")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-o", "--output", help="The output venn plot")
	parser.add_argument("-s", "--use_strand", action='store_true', default=False,
						help="Whether to use strand information when building keys")
	parser.add_argument("-m", "--multiclass", action='store_true', default=False,
						help="""Breakdown results into multiple classes:
						1) Matching intron
						2) Two matching splice sites but no matching intron (i.e. splice sites from different introns)
						3) One matching splice site
						4) No matching splice sites""")
	args = parser.parse_args()

	if not args.multiclass:

		ref_bed = Junction.createDict(args.reference, use_strand=args.use_strand)
		print("Loaded Reference file.  # junctions: ", len(ref_bed))

		# Load all bed files
		print("Results:")
		print("File\t#junc\t", Performance.shortHeader())

		recall = 0
		precision = 0
		f1 = 0

		for bf in args.input:
			bed_data = Junction.createDict(bf, use_strand=args.use_strand)

			# Build table
			tab = list()
			p = Performance()
			p.tp = len(ref_bed & bed_data)
			p.fp = len(bed_data - ref_bed)
			p.fn = len(ref_bed - bed_data)

			print(bf, "\t", len(bed_data), "\t", p)

			recall += p.recall()
			precision += p.precision()
			f1 += p.F1()

		if len(args.input) > 1:
			print("Mean recall: ", recall / len(args.input))
			print("Mean precision: ", precision / len(args.input))
			print("Mean f1: ", f1 / len(args.input))

		if not args.output == None and len(args.input) == 1:
			# Create Venns
			plt = figure(1, figsize=(6, 6))
			venn2(subsets=(p.fn, p.fp, p.tp), set_labels=(args.reference, args.input))
			plt.show()
			plt.savefig(args.output)

	else:

		ref_intron = Junction.createSet(args.reference, use_strand=args.use_strand)
		ref_ss = Junction.createSpliceSiteSet(args.reference, use_strand=args.use_strand)
		print("Loaded Reference BED file.  # introns: ", len(ref_intron), ", with ", len(ref_ss), " distinct splice sites")

		# Load all bed files
		print("Results:")
		print("Class 1 = Intron in ref")
		print("Class 2 = Both splice sites in ref")
		print("Class 3 = Only 1 splice site in ref")
		print("Class 4 = Novel")

		print("\t".join(["file", "class1", "class2", "class3", "class4"]))

		for jf in args.input:

			juncfileset = Junction.createSet(jf, use_strand=args.use_strand, keyonly=True)

			class1 = 0
			class2 = 0
			class3 = 0
			class4 = 0

			for j in juncfileset:
				key1 = j.startSplicsiteKey()
				key2 = j.endSplicesiteKey()

				if j in ref_intron:
					class1 += 1
				elif key1 in ref_ss and key2 in ref_ss:
					class2 += 1
				elif key1 in ref_ss or key2 in ref_ss:
					class3 += 1
				else:
					class4 += 1

			print("\t".join([jf, str(class1), str(class2), str(class3), str(class4)]))


main()
