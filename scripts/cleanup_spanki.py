#!/usr/bin/env python3

__author__ = 'maplesod'

import sys


class Junction:
	seq = ""
	start = 0
	end = 0
	strand = "."
	rest = ""

	def __init__(self):
		self.data = []

	def display(self):
		print(self.seq + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + self.strand + self.rest, end="")


junctions = list()

with open(sys.argv[1]) as f:
	# Skip header
	f.readline()

	for line in f:

		words = line.split("\t")

		parts1 = words[0].split(":")
		parts2 = parts1[1].split("_")

		rest = ""

		for i in range(1, len(words)):
			rest += "\t" + words[i]

		j = Junction()
		j.seq = parts1[0]
		j.start = int(parts2[0])
		j.end = int(parts2[1])
		j.strand = parts1[2]
		j.rest = rest
		junctions.append(j)

# Sort by
junctions.sort(key=lambda x: x.end)
junctions.sort(key=lambda x: x.start)
junctions.sort(key=lambda x: x.seq)

print(
	"refid\tstart\tend\tstrand\tdinucleotide\tintron_size\tannostatus\tgmcode\tregcode\tgeneassign\tgeneassignL\tgeneassignR\tunfilt_cov\tcov\tnormcov\toffsets\tentropy\thamming3\thamming5\tMAXmmes\tMAXminanc\tlirt\trirt\tirt\tdatrans\tdncov\tancov")

for j in junctions:
	j.display()
