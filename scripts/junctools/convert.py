import argparse
import sys
import collections

from .junction import JuncFactory, Junction, BedJunction


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


def loadgtf(filepath, dedup=False):
	'''
	Assumes that GTF is sorted (doesn't matter if it's sorted by transcript or exon first though)
	:param filepath:
	:return:
	'''


	# Create dict of transcripts with all exons in memory
	transcripts = collections.defaultdict(list)
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

	junctions = set() if dedup else []

	for t, exons in transcripts.items():
		le = None
		for i, e in enumerate(exons):
			if i > 0:
				j = BedJunction()
				j.refseq = e[0]
				j.start = int(le[2]) + 1
				j.end = int(e[1]) - 1
				j.strand = e[3]
				j.id = t

				if dedup:
					junctions.add(j)
				else:
					junctions.append(j)

			le = e

	if dedup:
		junctions = list(junctions)

	decstart(junctions)

	return junctions


def convert(args):
	in_type = JuncFactory[args.input_format.upper()]
	out_type = JuncFactory[args.output_format.upper()]

	# Check for any invalid input file formats
	if in_type == JuncFactory.GFF or in_type == JuncFactory.EGFF:
		raise "Currently we can only interpret GFF files containing intron features.  If you have this kind of file use \"igff\" as input type, or add introns to your gff."
	elif in_type == JuncFactory.EBED or in_type == JuncFactory.IBED or in_type == JuncFactory.TBED or in_type == JuncFactory.BED6:
		raise "Do not specify which type of BED file you wish to input, we will automatically detect and handle the BED file correctly.  Please rerun using \"bed\" as input type."

	# Check for any invalid output file formats
	if out_type == JuncFactory.GFF:
		raise "Must specify which particular type of GFF file to output: EGFF; IGFF"
	elif out_type == JuncFactory.BED:
		raise "Must specify which particular type of BED file to output: EBED, IBED, TBED, BED6"
	elif out_type == JuncFactory.SPANKI:
		raise "Haven't implemented SPANKI output yet"

	# Setup output handle
	o = sys.stdout
	if args.output != sys.stdout:
		o = open(args.output, mode='w')

	# Output header
	header = JuncFactory.create_from_enum(out_type).file_header()
	if header and header != "":
		print(header, file=o)

	# Whether or not to load the input file into memory or stream
	loadall = True if args.sort or in_type.isStreamable() or in_type == JuncFactory.GTF else False

	junction_set = set()
	junctions = []
	index = args.index_start

	if in_type == JuncFactory.GTF:
		junctions = loadgtf(args.input, dedup=args.dedup)
		args.sort = True	# Make sure we sort the output
	else:
		with open(args.input) as f:
			for line in f:
				j = JuncFactory.create_from_enum(in_type, use_strand=not args.ignore_strand).parse_line(line)

				if j:

					if args.dedup:
						if j.key not in junction_set:
							junction_set.add(j.key())
							if loadall:
								junctions.append(j)
					else:
						if loadall:
							junctions.append(j)

					if not loadall:
						c = JuncFactory.create_from_enum(out_type, use_strand=not args.ignore_strand, junc_to_copy=j)

						if args.reindex:
							c.id = args.prefix + str(index)
							index += 1

						if out_type.isBed() or out_type.isGFF:
							c.style = out_type

						print(c, file=o)

	if loadall:

		if args.sort:
			Junction.sort(junctions)
		if args.reindex:
			Junction.reindex(junctions, prefix=args.prefix, start=args.index_start)

		for j in junctions:

			# Convert
			c = JuncFactory.create_from_enum(out_type, use_strand=not args.ignore_strand, junc_to_copy=j)

			if out_type.isBed() or out_type.isGFF:
				c.style = out_type

			print(c, file=o)

	if args.output != sys.stdout:
		o.close()


def add_options(parser):
	parser.formatter_class = argparse.RawTextHelpFormatter

	parser.add_argument("-if", "--input_format", required=True, help='''The format of the input file to convert.''')

	parser.add_argument("-of", "--output_format", required=True, help='''The output format.''')

	parser.add_argument("-o", "--output", default=sys.stdout,
						help="Output to this file.  By default we print to stdout.")

	parser.add_argument("-is", "--ignore_strand", action='store_true', default=False,
						help="Whether or not to ignore strand when creating a key for the junction")

	parser.add_argument("-d", "--dedup", action='store_true', default=False,
						help="Whether or not to remove duplicate junctions")

	parser.add_argument("-s", "--sort", action='store_true', default=False,
						help="Whether or not to sort the junctions.  Note that sorting requires all junctions to be loaded into memory first.  This maybe an issue for very large input files.")

	parser.add_argument("-r", "--reindex", action='store_true', default=False,
						help="Whether or not to reindex the output.  The index is applied after prefix.")

	parser.add_argument("--index_start", type=int, default=0,
						help="The starting index to apply if the user requested reindexing")

	parser.add_argument("--prefix", default="junc_",
						help="The prefix to apply to junction ids if the user requested reindexing")

	parser.add_argument("--source", default="portcullis",
						help="Only relevant if output is GFF format, use this option to set the source column in the GFF")

	parser.add_argument("input", help="The input file to convert")
