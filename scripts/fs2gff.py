#!/usr/bin/env python
__author__ = 'maplesod'

import sys

with open(sys.argv[1]) as f:
    content = f.readline()

    index = 0;
    for line in f:
        words = line.split()

        print words[0] + "\tfinesplice\tjunction\t" + str(words[1]) + "\t" + str(words[2]) + "\t" + str(words[3]) + "\t.\t.\tID=junc_" + str(index) + ";"
        index += 1