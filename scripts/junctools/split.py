#!/usr/bin/env python3
__author__ = 'maplesod'

from .junction import *
from .performance import Performance


def write(filepath, input, juncset):
	with open(filepath, "wt") as out:
		header = JuncFactory.create_from_file(filepath).file_header()
		print(header, file=out)
		with open(input) as fin:
			for line in fin:
				junc = JuncFactory.create_from_file(input).parse_line(line, fullparse=False)
				if junc:
					if junc.key in juncset:
						oj = JuncFactory.create_from_file(filepath, junc_to_copy=junc)
						print(oj, file=out)


def split(args):
	print("\t".join(["File", "distinct", "total"]))
	# Load reference bed file as a set
	ref_juncs, ref_count = Junction.createJuncSet(args.reference)
	print("\t".join(["reference", str(len(ref_juncs)), str(ref_count)]))
	pass_juncs, pass_count = Junction.createJuncSet(args.passfile)
	print("\t".join(["pass", str(len(pass_juncs)), str(pass_count)]))
	fail_juncs, fail_count = Junction.createJuncSet(args.failfile)
	print("\t".join(["fail", str(len(fail_juncs)), str(fail_count)]))
	print()

	filename, ext = os.path.splitext(args.passfile)
	filename, ext2 = os.path.splitext(args.failfile)

	if ext != ext2:
		raise ValueError("Pass and fail files should be the same format")

	tp = pass_juncs & ref_juncs
	tn = fail_juncs - ref_juncs
	fp = pass_juncs - ref_juncs
	fn = ref_juncs - fail_juncs

	p = Performance(tp=len(tp), fp=len(fp), tn=len(tn), fn=len(fn))
	print(p.longHeader())
	print(p.longStr())

	# Filenames
	tpfile = args.output_prefix + ".TP" + ext
	fnfile = args.output_prefix + ".FN" + ext
	fpfile = args.output_prefix + ".FP" + ext
	tnfile = args.output_prefix + ".TN" + ext

	write(tpfile, args.passfile, tp)
	write(tnfile, args.failfile, tn)
	write(fpfile, args.passfile, fp)
	write(fnfile, args.reference, fn)

	print("Output files saved to:", args.output_prefix)


def add_options(parser):
	parser.add_argument("reference", help="The reference junction file")
	parser.add_argument("passfile", help="The junction file containing junctions that pass a filter")
	parser.add_argument("failfile", help="The junction file containing junctions failing a filter")
	parser.add_argument("-o", "--output_prefix", default="split",
						help="Prefix for output files")
