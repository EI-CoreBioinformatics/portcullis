#!/usr/bin/env python3
__author__ = 'maplesod'

import sys


class Junction:
	seq = ""
	start = 0
	end = 0
	strand = "."
	cov = 0

	def __init__(self):
		self.data = []

	def display(self):
		print(self.seq + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + self.strand + "\t" + str(self.cov))

	def key(self):
		return (self.seq, self.start, self.end, self.strand)


with open(sys.argv[1]) as f:
	junctions = [];
	junction_set = {}

	curr_transcript_seq = ""
	curr_transcript_start = 0
	curr_transcript_end = 0
	curr_transcript_strand = ""

	last_exon_seq = ""
	last_exon_start = 0
	last_exon_end = 0
	last_exon_strand = ""

	index = 0;
	for line in f:

		if not line.startswith("#"):
			words = line.split()

			start = int(words[3])
			end = int(words[4])

			dif = end - start

			if words[2] == "transcript":

				curr_transcript_seq = words[0]
				curr_transcript_start = start
				curr_transcript_end = end
				curr_transcript_strand = words[6]

				# Wipe last exon
				last_exon_seq = ""
				last_exon_start = 0
				last_exon_end = 0
				last_exon_strand = ""

			elif words[2] == "exon":

				# Double check we are still in the current transcript and there has already been an exon in this transcript
				if curr_transcript_end >= last_exon_end and curr_transcript_seq == last_exon_seq and curr_transcript_strand == last_exon_strand:

					j = Junction()
					j.seq = words[0]
					j.start = last_exon_end + 1
					j.end = start - 1
					j.strand = words[6]

					if j.key() not in junction_set:
						junctions.append(j)
						junction_set[j.key()] = j

				last_exon_seq = words[0]
				last_exon_start = start
				last_exon_end = end
				last_exon_strand = words[6]

# Sort by
junctions.sort(key=lambda x: x.end)
junctions.sort(key=lambda x: x.start)
junctions.sort(key=lambda x: x.seq)

print("track name=\"junctions\"")

index = 0;
for j in junctions:
	# Print bed line adjusting offsets as necessary to do the correct conversion
	print(j.seq + "\t" + str(j.start - 2) + "\t" + str(j.end + 1) + "\tjunc_" + str(index) + "\t" + str(
		j.cov) + "\t" + j.strand + "\t" + str(j.start - 1) + "\t" + str(j.end) + "\t255,0,0\t2\t1,1\t0," + str(
		j.end - j.start))

	index += 1
