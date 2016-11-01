from matplotlib_venn import venn2
from pylab import figure

from .performance import Performance
from .junction import *


def compare(args):
	if not args.multiclass:

		rp = set()
		rn = set()
		ref_entries = 0

		if args.labels:
			rp, ref_entries = Junction.createJuncSet(args.reference[0], use_strand=args.use_strand)
		else:
			rp, rn, ref_entries = Junction.createMarkedupJuncSets(args.reference[0], use_strand=args.use_strand)

		print()
		print("Reference:")
		print(" - # total junctions:", ref_entries)
		if args.labels:
			print(" - # distinct positive junctions:", len(rp))
			print(" - # distinct negative junctions:", len(rn))
		else:
			print(" - # distinct junctions:", len(rp))

		print()

		# Load all junction files
		print("\t".join(
			["File", "distinct", "total", Performance.longHeader() if args.labels else Performance.shortHeader()]))

		recall = 0
		precision = 0
		f1 = 0

		for f in args.input:
			junc_set, bed_entries = Junction.createJuncSet(f, use_strand=args.use_strand)

			# Build table
			tab = list()
			p = Performance()

			if args.labels:
				p.tp = len(junc_set & rp)
				p.fp = len(junc_set & rn)
				p.fn = len(rp - junc_set)
				p.tn = len(rn - junc_set)
			else:
				p.tp = len(rp & junc_set)
				p.fp = len(junc_set - rp)
				p.fn = len(rp - junc_set)

			print("\t".join([f, str(len(junc_set)), str(bed_entries), p.longStr() if args.labels else str(p)]))

			recall += p.recall()
			precision += p.precision()
			f1 += p.F1()

		if len(args.input) > 1:
			print("Mean recall: ", recall / len(args.input))
			print("Mean precision: ", precision / len(args.input))
			print("Mean f1: ", f1 / len(args.input))

	else:

		ref_set, ref_entries = Junction.createJuncSet(args.reference[0], use_strand=args.use_strand)
		ref_ss = Junction.createSpliceSiteSet(args.reference[0], use_strand=args.use_strand)
		print()
		print("Reference:")
		print(" - # distinct junctions:", len(ref_set))
		print(" - # total junctions:", ref_entries)
		print(" - # distinct splice sites:", len(ref_ss))
		print()

		# Load all bed files
		print("Result legend:")
		print("Class 1 = Intron in ref")
		print("Class 2 = Both splice sites in ref")
		print("Class 3 = Only 1 splice site in ref")
		print("Class 4 = Novel")
		print()

		print("\t".join(["file", "class1", "class2", "class3", "class4"]))

		for jf in args.input:

			juncs, entries = Junction.createDict(jf, use_strand=args.use_strand, fullparse=True)

			class1 = 0
			class2 = 0
			class3 = 0
			class4 = 0

			for key, value in juncs.items():
				key1 = value.startSplicesiteKey()
				key2 = value.endSplicesiteKey()

				if value.key in ref_set:
					class1 += 1
				elif key1 in ref_ss and key2 in ref_ss:
					class2 += 1
				elif key1 in ref_ss or key2 in ref_ss:
					class3 += 1
				else:
					class4 += 1

			print("\t".join([jf, str(class1), str(class2), str(class3), str(class4)]))


def add_options(parser):
	parser.add_argument("reference", nargs=1, help="The junction file to treat as the reference")
	parser.add_argument("input", nargs="+", help="One or more junction files to compare against the reference")
	parser.add_argument("-s", "--use_strand", action='store_true', default=False,
						help="Whether to use strand information when building keys")
	parser.add_argument("-l", "--labels",
						help="Path to a file containing labels for the reference indicating whether or not each reference junction is genuine (as generated using the markup tool).  If provided this script produces a much richer performance analysis.  Not compatible with '--multiclass'")
	parser.add_argument("-m", "--multiclass", action='store_true', default=False,
						help="""Breakdown results into multiple classes:
						1) Matching intron
						2) Two matching splice sites but no matching intron (i.e. splice sites from different introns)
						3) One matching splice site
						4) No matching splice sites""")
