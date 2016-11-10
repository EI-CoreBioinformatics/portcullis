#!/usr/bin/env python3

from .junction import *


def markup(args):
	ref_set, ref_entries = Junction.createJuncSet(args.reference[0], use_strand=args.use_strand)
	print()
	print("Reference:")
	print(" - # distinct junctions:", len(ref_set))
	print(" - # total junctions:", ref_entries)
	print()
	print("\t".join(["File", "Total", "InRef", "OutRef"]))

	nb_pos = 0
	nb_neg = 0
	for filepath in args.input:

		head, tail = os.path.split(filepath)
		outfile = (filepath + ".res") if not args.output_dir else (os.path.join(args.output_dir, tail + ".res"))

		res = open(outfile, "w")

		with open(filepath) as f:
			for line in f:
				junc = JuncFactory.create_from_file(filepath, use_strand=args.use_strand).parse_line(line,
																									 fullparse=False)
				if junc:

					if junc.key in ref_set:
						print("1", file=res)
						nb_pos += 1
					else:
						print("0", file=res)
						nb_neg += 1

		res.close()

		print("\t".join([filepath, str(nb_pos + nb_neg), str(nb_pos), str(nb_neg)]))


def add_options(parser):
	parser.add_argument("reference", nargs=1, help="The junction file to treat as the reference")
	parser.add_argument("input", nargs="+", help="One or more junction files to compare against the reference")
	parser.add_argument("-o", "--output_dir", help='''If output dir is specified this will create output files for each
					input file with a .res extension indicating whether or not the junction was found in the reference.
					By default we write out a .res file in the same directory as the input file was found in.''')
	parser.add_argument("-s", "--use_strand", action='store_true', default=False,
						help="Whether to use strand information when building keys")
