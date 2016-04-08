#!/usr/bin/env python3

import itertools
import argparse
import bed12
import performance
from matplotlib_venn import venn2
from pylab import figure


def main():
	parser = argparse.ArgumentParser("Script to compare bed file against reference bed")
	parser.add_argument("input", help="The BED file to analyse")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-o", "--output", help="The output venn plot")
	args = parser.parse_args()

	ref_bed = bed12.loadbed(args.reference, False, False)
	print("Loaded Reference BED file.  # junctions: ", len(ref_bed))

	# Load all bed files
	bed_data = bed12.loadbed(args.input, False, False)
	print("Loaded: " + args.input + "; # junctions: ", len(bed_data))

	# Build table
	tab = list()
	p = performance.PEntry()
	p.tp = len(ref_bed & bed_data)
	p.fp = len(bed_data - ref_bed)
	p.fn = len(ref_bed - bed_data)

	print("Results:")
	print(performance.PEntry.header())
	print(p)

	if not args.output == None:
		# Create Venns
		plt = figure(1, figsize=(6, 6))
		venn2(subsets=(p.fn, p.fp, p.tp), set_labels=(args.reference, args.input))
		plt.show()
		plt.savefig(args.output)


main()
