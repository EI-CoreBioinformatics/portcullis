#!/usr/bin/env python3

import sys
import os

with open(sys.argv[1]) as f:

    # Skip header
    header = f.readline()
    print(header, end="")

    for line in f:

        words = line.split("\t")

        entropy = float(words[13])
        hamming3 = int(words[14])
        hamming5 = int(words[15])

        # These are the filters described in the SPANKI paper, which don't involve using the reference
        if entropy > 2.0 and hamming3 >= 8 and hamming5 >= 8:
            print(line, end="")