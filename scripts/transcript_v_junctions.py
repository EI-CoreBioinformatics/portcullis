#!/usr/bin/env python3
__author__ = 'maplesod'

import sys
import argparse
import os
import bed12


def main():
	parser = argparse.ArgumentParser("Script to determine how many transcripts have correct intron chains")
	parser.add_argument("transcripts_in", help="The GTF file containing transcripts")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	parser.add_argument("-b", "--bed", help="If provided this is a representation of valid junctions in the provided reference which can be used for filtering")
	#parser.add_argument("-o", "--output", required=True, help="The directory")
	args = parser.parse_args()

	#if not os.path.exists(args.output):
	#	os.mkdir(args.output)

	#os.system("gtf2bed " + args.transcript + " > " + args.output + "/input.bed")
	#print("Extracted junctions from transcripts")

	input_bed = {}

	if args.bed:
		input_bed = bed12.loadbed(args.bed, False, False)
		print("Loaded input BED file.  # junctions: ", len(input_bed))

	ref_bed = bed12.loadbed(args.reference, False, False)
	print("Loaded Reference BED file.  # junctions: ", len(ref_bed))


	with open(args.transcripts_in, 'r') as f:
		junctions = []
		junction_set = {}

		curr_transcript_seq = ""
		curr_transcript_start = 0
		curr_transcript_end = 0
		curr_transcript_strand = ""

		last_exon_seq = ""
		last_exon_start = 0
		last_exon_end = 0
		last_exon_strand = ""

		total_transcripts = 0;
		valid_transcripts = 0
		invalid_transcripts = 0
		valid_mono_exonic_transcripts = 0
		ignored_transcripts = 0

		total_junctions = 0;
		valid_junctions = 0
		invalid_junctions = 0


		for line in f:

			if not line.startswith("#"):
				words = line.split()

				start = int(words[3])
				end = int(words[4])

				dif = end - start

				if words[2] == "transcript":
					total_transcripts+=1

					# If this is a new transcript then check if we have some introns to work with.  If so then validate
					if len(junctions) > 0:
						valid = True
						ignore = False
						for j in junctions:
							if args.bed and j not in input_bed:
								ignore = True
								break
							elif j not in ref_bed:
								valid = False
								break
						if valid and not ignore:
							valid_transcripts += 1
						elif ignore:
							ignored_transcripts += 1
						else:
							invalid_transcripts += 1
					elif len(junctions) == 0:
						valid_mono_exonic_transcripts += 1

					curr_transcript_seq = words[0]
					curr_transcript_start = start
					curr_transcript_end = end
					curr_transcript_strand = words[6]

					# Wipe last exon
					last_exon_seq = ""
					last_exon_start = 0
					last_exon_end = 0
					last_exon_strand = ""

					# Wipe intron chain
					junctions.clear()

				elif words[2] == "exon":

					# Double check we are still in the current transcript and there has already been an exon in this transcript
					if curr_transcript_end >= last_exon_end and curr_transcript_seq == last_exon_seq and curr_transcript_strand == last_exon_strand:

						j = bed12.BedEntry(use_strand=False)
						j.chrom = words[0]
						j.thick_start = last_exon_end
						j.thick_end = start - 1
						j.strand = last_exon_strand
						#print(j)

						total_junctions += 1
						junctions.append(j)
						if j in ref_bed:
							valid_junctions += 1
						else:
							invalid_junctions += 1

						if j not in junction_set:
							junction_set[j] = j

					last_exon_seq = words[0]
					last_exon_start = start
					last_exon_end = end
					last_exon_strand = words[6]

		if len(junctions) > 0:
			for j in junctions:
				valid = True
				if j not in ref_bed:
					valid = False
					break
				if valid:
					valid_transcripts += 1
				else:
					invalid_transcripts += 1
		elif len(junctions) == 0:
			valid_transcripts += 1

		print("Total transcripts", total_transcripts)
		print("Valid transcripts", valid_transcripts)
		print("Invalid transcripts", invalid_transcripts)
		print("Mono-exonic transcripts", valid_mono_exonic_transcripts)
		print("Ignored transcripts", ignored_transcripts)
		print("Total valid", valid_transcripts + valid_mono_exonic_transcripts)


		print("Total junctions", total_junctions)
		print("Valid junctions", valid_junctions)
		print("Invalid junctions", invalid_junctions)

main()