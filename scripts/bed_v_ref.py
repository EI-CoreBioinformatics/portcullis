#!/usr/bin/env python3

import itertools
import argparse
import bed12
from performance import Performance
from matplotlib_venn import venn2
from pylab import figure


def main():
	parser = argparse.ArgumentParser("Script to compare bed file against reference bed")
	parser.add_argument("input", nargs="+", help="The BED file to analyse")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-o", "--output", help="The output venn plot")
	args = parser.parse_args()

	ref_bed = bed12.loadbed(args.reference, False, False)
	print("Loaded Reference BED file.  # junctions: ", len(ref_bed))

	# Load all bed files
	print("Results:")
	print("File\t#junc\t", Performance.shortHeader())

	recall = 0
	precision = 0
	f1 = 0

	for bf in args.input:
		bed_data = bed12.loadbed(bf, False, False)

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


main()
