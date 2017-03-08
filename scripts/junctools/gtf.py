import argparse
import collections
from enum import Enum, unique

from .junction import Junction, BedJunction

@unique
class Mode(Enum):
	MARKUP = 1
	FILTER = 2

def decstart(junctions):
	last_id = ""
	index = 1
	for j in junctions:
		j.start -= 1
		j.end -= 1

		if last_id == "" or last_id != j.id:
			last_id = j.id
			index = 1

		j.id += "_junc" + str(index)
		index += 1

def loadgtf(filepath, use_strand=False):
	'''
	Assumes that GTF is sorted (doesn't matter if it's sorted by transcript or exon first though)
	:param filepath:
	:return:
	'''


	# Create dict of transcripts with all exons in memory
	transcripts = collections.defaultdict(list)
	nb_transcripts = 0
	with open(filepath) as f:
		for line in f:
			if not line.startswith("#"):
				parts = line.split('\t')
				if len(parts) == 9 and parts[2] == "exon":
					tags = parts[8].split(';')
					for tag in tags:
						t = tag.strip()
						if t:
							tag_parts = t.split()
							if tag_parts[0] == "transcript_id":
								t = tag_parts[1].strip()
								transcript_id = t[1:-1] if t[0] == '\"' else t
								transcripts[transcript_id].append([parts[0], parts[3], parts[4], parts[6]])
				elif len(parts) == 9 and parts[2] == "transcript":
					nb_transcripts += 1

	intron_chains = collections.defaultdict(list)
	junc_set = set()
	nb_introns = 0

	for t, exons in transcripts.items():
		le = None
		for i, e in enumerate(exons):
			if i > 0:
				j = BedJunction(use_strand=use_strand)
				j.refseq = e[0]
				j.start = int(le[2]) + 1
				j.end = int(e[1]) - 1
				j.strand = e[3]
				j.id = t

				intron_chains[t].append(j)
				junc_set.add(j)
				nb_introns += 1

			le = e

	last_id = ""
	index = 1
	for t, introns in intron_chains.items():
		for j in introns:
			j.start -= 1
			j.end -= 1

			if last_id == "" or last_id != j.id:
				last_id = j.id
				index = 1

			j.id += "_junc" + str(index)
			index += 1

	#for j in junc_set:
	#	j.start -= 1
	#	j.end -= 1

	junc_key_set = set()
	for j in junc_set:
		junc_key_set.add(j.key)

	return intron_chains, junc_key_set, nb_transcripts, nb_introns

def gtf(args):
	mode = Mode[args.mode.upper()]
	print("Running junctools GTF in", mode.name, "mode")

	print("Loading junctions ...",end="")
	port_juncs, port_count = Junction.createJuncSet(args.junctions, use_strand=not args.ignore_strand)
	print(" done.  Found", port_count, "junctions.")
	print("Loading transcripts ...",end="")
	intron_chains, junc_set, nb_transcripts, nb_introns = loadgtf(args.input, use_strand=not args.ignore_strand)
	print(" done.")
	nb_monoexonic = nb_transcripts - len(intron_chains)
	nb_multiexonic = len(intron_chains)
	print("Found", nb_transcripts, "transcripts. ", nb_monoexonic, "are monoexonic and", nb_multiexonic, "are multiexonic.")
	print("Found", len(junc_set), "distinct junctions and", nb_introns, "total junctions in GTF.")
	print()

	print("Doing junction level comparison ...", end="")
	juncs_inport = port_juncs & junc_set
	print(" done.")
	recall = (len(juncs_inport) / port_count) * 100.0
	precision = (len(juncs_inport) / len(junc_set)) * 100.0
	print(len(juncs_inport), "/", port_count, "(" + "{0:.2f}".format(recall) + ") of valid junctions.")
	print(len(juncs_inport), "/", len(junc_set), "(" + "{0:.2f}".format(precision) + ") junctions supported by input.")
	print()


	print("Identifying invalid multi-exonic transcripts ... ", end="")
	invalid_transcripts = collections.defaultdict(list)
	for t, introns in intron_chains.items():
		for index, j in enumerate(introns):
			if not j.key in port_juncs:
				invalid_transcripts[t].append(str(j.start+1) + "_" + str(j.end+1))
	print(" done.  Found", len(invalid_transcripts), "/", nb_multiexonic, "invalid multi-exonic transcripts.")
	print()


	print("Writing output to", args.output," ... ", end="")
	o = open(args.output, mode='w')

	# Create dict of transcripts with all exons in memory
	with open(args.input) as f:
		for l in f:
			line = l.strip()
			if not line.startswith("#"):
				parts = line.split('\t')
				if len(parts) == 9 and (parts[2] == "exon" or parts[2] == "transcript"):
					tags = parts[8].split(';')
					for tag in tags:
						t = tag.strip()
						if t:
							tag_parts = t.split()
							if tag_parts[0] == "transcript_id":
								tid = tag_parts[1].strip()
								transcript_id = tid[1:-1] if tid[0] == '\"' else tid
								if transcript_id in invalid_transcripts:
									if mode != Mode.FILTER:
										if mode == Mode.MARKUP and parts[2] == "transcript":
											bad_juncs = ",".join(invalid_transcripts[transcript_id])
											o.write(line + " introns \"invalid(" + bad_juncs + ")\";\n")
										else:
											o.write(line + "\n")
								else:
									if mode == Mode.MARKUP and parts[2] == "transcript":
										o.write(line + " introns \"valid\";\n")
									else:
										o.write(line + "\n")
				else:
					o.write(line + "\n")
			else:
				o.write(line + "\n")

	o.close()
	print(" done.")


def add_options(parser):
	parser.formatter_class = argparse.RawTextHelpFormatter

	parser.add_argument("-is", "--ignore_strand", action='store_true', default=False,
						help="Whether or not to ignore strand when looking for junctions")

	parser.add_argument("-j", "--junctions", required=True, help='''The file containing junctions that should be found in the GTF.''')

	parser.add_argument("-o", "--output", required=True, default="junctools.out.gtf",
						help="The filtered or markedup GTF output.  By default we print to stdout.")

	parser.add_argument("mode", help='''GTF operation to apply.  See above for details.  Available options:
 - filter
 - markup
''')

	parser.add_argument("input", help="The input GTF file to convert")
