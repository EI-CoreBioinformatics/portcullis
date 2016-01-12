#!/usr/bin/env python3

# Script is intended to be run on cufflinks "isoform.tracking" output file when cufflinks is run using the --GTF option.  Output is a two column tab delimited file containing the transcript id in the first column and integer rounded coverage in the second column.  This file can then be passed to spanki_transcripts using the "-t" option, producing a set of simulated transcripts based on the reference transcript set but at levels of coverage.

import sys
import os

with open(sys.argv[1]) as f:

    # Skip header
    f.readline()
    print ("transcript_id\tcov")

    for line in f:

        words = line.split("\t")

        id = words[0]
        cov = int(float(words[8]))

        print (id + "\t" + str(cov))
