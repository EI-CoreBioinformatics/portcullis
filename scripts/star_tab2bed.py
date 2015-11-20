#!/usr/bin/env python3

import sys
import os


class Junction:

    seq = ""
    start = 0
    end = 0
    strand = "."
    cov = 0
    
    def __init__(self):
        self.data = []

    def display(self):
        print (self.seq + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + self.strand + "\t" + str(self.cov))



junctions = list()

with open(sys.argv[1]) as f:

    # Skip header
    f.readline()

    for line in f:


        words = line.split("\t")

        j = Junction()
        j.seq = words[0]
        j.start = int(words[1])
        j.end = int(words[2])
        j.strand = "+" if words[3] == "1" else "-" if words[3] == "2" else "?" 
        j.cov = int(words[6])
        junctions.append(j)

# Sort by
junctions.sort(key=lambda x: x.end)
junctions.sort(key=lambda x: x.start)
junctions.sort(key=lambda x: x.seq)


#for j in junctions:
#    j.display()

print ("track name=\"junctions\"")

index = 0;
for j in junctions:
    print (j.seq + "\t" + str(j.start - 2) + "\t" + str(j.end + 1) + "\tjunc_" + str(index) + "\t" + str(j.cov) + "\t" + j.strand + "\t" + str(j.start - 1) + "\t" + str(j.end) + "\t255,0,0\t2\t1,1\t0," + str(j.end - j.start))

    index += 1

