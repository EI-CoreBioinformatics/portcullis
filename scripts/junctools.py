#!/usr/bin/env python3

import argparse
import sys

import junctools.compare
import junctools.convert
import junctools.markup
import junctools.set

def main():
	call_args = sys.argv[1:]

	parser = argparse.ArgumentParser(
		"""This script contains a number of tools for manipulating junction files.""",
		formatter_class=argparse.RawTextHelpFormatter)
	subparsers = parser.add_subparsers(
		title="Junction tools")

	compare_parser = subparsers.add_parser("compare",
										   help="Compares junction files.")
	junctools.compare.add_options(compare_parser)
	compare_parser.set_defaults(func=junctools.compare.compare)

	convert_parser = subparsers.add_parser("convert",
										   help="Converts junction files between various formats.")
	junctools.convert.add_options(convert_parser)
	convert_parser.set_defaults(func=junctools.convert.convert)

	markup_parser = subparsers.add_parser("markup",
										 help="Marks whether each junction in the input can be found in the reference or not.")
	junctools.markup.add_options(markup_parser)
	markup_parser.set_defaults(func=junctools.markup.markup)

	merge_parser = subparsers.add_parser("set",
										 help="Apply set operations to two or more junction files.")
	junctools.set.add_options(merge_parser)
	merge_parser.set_defaults(func=junctools.set.setops)

	args = parser.parse_args(call_args)
	if hasattr(args, "func"):
		args.func(args)
	else:
		parser.print_help()


if __name__ == '__main__':
	main()
