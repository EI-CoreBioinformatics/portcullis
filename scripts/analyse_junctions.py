#!/usr/bin/env python3

import os
import argparse
import bed12
from performance import Performance


def main():
	parser = argparse.ArgumentParser("Script to create the Venn Plots from BED files")
	parser.add_argument("input", nargs='+', help="The directory containing BED files from pipeline")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-o", "--output", required=True, help="The output prefix")
	parser.add_argument("-f", "--filter", )
	args = parser.parse_args()

	ref_bed = bed12.loadbed(args.reference, False, False)
	print("Loaded Reference BED file.  # junctions: " + str(len(ref_bed)))

	# Load all bed files
	bed_data = {}
	aligners = set()
	reads = set()
	junc_analysers = set()
	for bed_file in args.input:
		bed_base = os.path.splitext(os.path.basename(bed_file))[0]
		parts = bed_base.split('-')
		if (not parts[0] == "trinity"):
			aligners.add(parts[0])
			reads.add(parts[1])
			junc_analysers.add(parts[2])
			bed_data[parts[2] + "-" + parts[0]] = bed12.loadbed(bed_file, False, False)
			print("Loaded: " + bed_file + "; # junctions: " + str(len(bed_data[parts[2] + "-" + parts[0]])))

	print("Found these aligners: " + ', '.join(aligners))
	print("Found these reads: " + ', '.join(reads))
	print("Found these junction analysis tools: " + ', '.join(junc_analysers))

	# Build table
	tab = []
	for a in aligners:
		for j in junc_analysers:
			p = Performance()
			p.tp = len(bed_data[j + "-" + a] & ref_bed)
			p.fp = len(bed_data[j + "-" + a] - ref_bed)
			p.fn = len(ref_bed - bed_data[j + "-" + a])
			tab.append(a + "\t" + j + "\t" + p.__str__())

	# Output table to disk
	with open(args.output + "-junc_analysis.tab", "w") as tab_out:
		print("Aligner\tFilter\t" + Performance.shortHeader(), file=tab_out)
		for p in tab:
			print(p, file=tab_out)


main()
