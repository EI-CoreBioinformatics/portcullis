#!/usr/bin/env python3

import sys
import bed12

with open(sys.argv[1]) as f:
	# Skip header
	f.readline()
	print("track name=\"junctions\"")

	for line in f:
		print(bed12.BedEntry.create_from_line(line, use_strand=True, tophat=True))
