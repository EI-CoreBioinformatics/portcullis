#!/usr/bin/env python3

import sys
import os
import argparse


def main():
	parser = argparse.ArgumentParser("Script to combine tab files together")
	parser.add_argument("input", nargs="+", help="A list of tab files to merge")
	args = parser.parse_args()

	# Load tab file and produce matrix
	tab = []
	count = 0

	# Print header
	with open(args.input[0]) as f:
		print(f.readline(), end="")

	for i in args.input:
		with open(i) as f:
			# Skip header
			f.readline()

			for line in f:
				line.strip()
				if len(line) > 1:
					print(line, end="")
					count += 1

	print("# tab entries: " + str(count) + " from " + str(len(args.input)) + " input files", file=sys.stderr)


main()
