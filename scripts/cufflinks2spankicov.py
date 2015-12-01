#!/usr/bin/env python3

import sys
import os

with open(sys.argv[1]) as f:

    # Skip header
    f.readline()
    print ("transcript_id\tcov")

    for line in f:

        words = line.split("\t")

        id = words[0]
        cov = float(words[8])

        print (id + "\t" + str(cov))
