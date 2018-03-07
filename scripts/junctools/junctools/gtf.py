import argparse
import collections
import sys
from enum import Enum, unique

from .junction import Junction, BedJunction
from .performance import Performance

@unique
class Mode(Enum):
	MARKUP = 1
	FILTER = 2
	COMPARE = 3

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

	intron_chains = collections.defaultdict(list)
	junc_set = set()
	nb_introns = 0
	monoexonics = set()

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
		if len(exons) == 1:
			# HACK! Use a junction to represent the mono-exonic transcript
			e = exons[0]
			j = BedJunction(use_strand=use_strand)
			j.refseq = e[0]
			j.start = int(e[1])
			j.end = int(e[2])
			j.strand = e[3]
			j.id = t
			monoexonics.add(j)


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

	return intron_chains, junc_key_set, len(transcripts), nb_introns, monoexonics

def run_compare(args, ref_juncs, ref_monos, ref_ics):
	print("\t".join(["file", "j_distinct", "j_total", "j_tp", "j_fp", "j_fn", "j_recall", "j_precision", "j_f1",
					 "t_transcripts", "t_monoexonic", "t_multiexonic", "t_supported", "t_unsupported", "t_precision",
					 "mt_tp", "mt_fp", "mt_fn", "mt_recall", "mt_precision", "mt_f1"
					 "ic_tp", "ic_fp", "ic_fn", "ic_recall", "ic_precision", "ic_f1"]))
	sys.stdout.flush()
	for i in args.input:
		intron_chains, junc_set, nb_transcripts, nb_introns, monoexonics = loadgtf(i, use_strand=not args.ignore_strand)
		nb_monoexonic = nb_transcripts - len(intron_chains)

		if nb_monoexonic != len(monoexonics):
			print("Monoexonic numbers don't agree for", i, nb_monoexonic, len(monoexonics), file=sys.stderr)

		nb_multiexonic = len(intron_chains)
		if nb_multiexonic <= 0:
			print("skipped...nb_multiexonic=0")
			continue

		jr_tp = len(ref_juncs & junc_set)
		jr_fn = len(ref_juncs - junc_set)
		jr_fp = len(junc_set - ref_juncs)
		jr_perf = Performance(tp=jr_tp, fn=jr_fn, fp=jr_fp)

		nb_in_ref = 0
		ic_tp = 0
		ic_fp = 0

		for t, introns in intron_chains.items():

			is_valid_by_ref = True

			# Start key to be used to represent this intron chain
			key = introns[0].refseq + "_" + introns[0].strand

			# Loop through introns in this transcript
			for index, j in enumerate(introns):
				in_ref = j.key in ref_juncs

				# Add intron start and stop to intron chain key
				key += "_" + str(j.start) + "_" + str(j.end)

				if not in_ref:
					is_valid_by_ref = False
					break

			if is_valid_by_ref:
				nb_in_ref += 1

			# Check if intron chain is found
			if key in ref_ics:
				ic_tp += 1
			else:
				ic_fp += 1

		mt_tp = len(ref_monos & monoexonics)
		mt_fn = len(ref_monos - monoexonics)
		mt_fp = len(monoexonics - ref_monos)
		mt_perf = Performance(tp=mt_tp, fp=mt_fp, fn=mt_fn)

		nb_unsupported = nb_multiexonic - nb_in_ref
		t_precision = (nb_in_ref / nb_multiexonic) * 100.0

		ic_fn = len(ref_ics) - ic_tp
		ic_perf = Performance(tp=ic_tp, fn=ic_fn, fp=ic_fp)

		print("\t".join(str(_) for _ in [i, len(junc_set), nb_introns,
										 jr_tp, jr_fp, jr_fn,
										 "{0:.2f}".format(jr_perf.recall()), "{0:.2f}".format(jr_perf.precision()), "{0:.2f}".format(jr_perf.F1()),
										 nb_transcripts, nb_monoexonic, nb_multiexonic, nb_in_ref, nb_unsupported, "{0:.2f}".format(t_precision),
										 mt_tp, mt_fp, mt_fn,
										 "{0:.2f}".format(mt_perf.recall()), "{0:.2f}".format(mt_perf.precision()), "{0:.2f}".format(mt_perf.F1()),
										 ic_tp, ic_fp, ic_fn,
										 "{0:.2f}".format(ic_perf.recall()), "{0:.2f}".format(ic_perf.precision()), "{0:.2f}".format(ic_perf.F1())]))
		sys.stdout.flush()

def keyFromIC(ic):
	if len(ic) > 0:
		key = ic[0].refseq + "_" + ic[0].strand
		for i in ic:
			key += "_" + str(i.start) + "_" + str(i.end)
		return key
	else:
		return None

def convert_ic_map(ic):
	mod = set()
	for t, introns in ic.items():
		key = keyFromIC(introns)
		if key:
			mod.add(key)

	return mod



def gtf(args):
	mode = Mode[args.mode.upper()]
	print("# Running junctools GTF in", mode.name, "mode")

	if mode != Mode.COMPARE and (len(args.input) > 1 or len(args.input) == 0):
		raise SyntaxError("This mode can only take a single GTF file as input")

	if args.junctions == None and args.transcripts == None:
		raise SyntaxError("You must either specify a file containing junctions (-j) or transcripts (-t) to use as a reference against your input files(s).")
	elif args.junctions and args.transcripts:
		raise SyntaxError(
			"You must either specify a file containing junctions (-j) or transcripts (-t) to use as a reference against your input files(s). Not both.")


	if args.junctions:
		print("# Loading junctions from reference junctions file ...",end="")
		ref_juncs, ref_junc_count = Junction.createJuncSet(args.junctions, use_strand=not args.ignore_strand)
		print(" done.  Found", ref_junc_count, "junctions (", len(ref_juncs), "distinct ) .")
	else:
		print("# Extracting junctions from reference transcript file ...", end="")
		ref_intron_chains, ref_juncs, ref_transcript_count, ref_junc_count, ref_monos = loadgtf(args.transcripts, use_strand=not args.ignore_strand)
		ref_ics = convert_ic_map(ref_intron_chains)
		print(" done.  Found", ref_junc_count, "junctions (", len(ref_juncs), "distinct ) and", ref_transcript_count, "transcripts (", len(ref_ics), "distinct intron chains).")

	if mode == Mode.COMPARE:
		run_compare(args, ref_juncs, ref_monos, ref_ics)

	else:
		print("Loading transcripts ...",end="")
		intron_chains, junc_set, nb_transcripts, nb_introns, monoexonics = loadgtf(args.input[0], use_strand=not args.ignore_strand)
		print(" done.")
		nb_monoexonic = nb_transcripts - len(intron_chains)
		nb_multiexonic = len(intron_chains)
		print("Found", nb_transcripts, "transcripts. ", nb_monoexonic, "are monoexonic and", nb_multiexonic, "are multiexonic.")
		print("Found", len(junc_set), "distinct junctions and", nb_introns, "total junctions in GTF.")
		print()

		print("Doing junction level comparison ...", end="")
		juncs_inport = ref_juncs & junc_set
		print(" done.")
		recall = (len(juncs_inport) / ref_junc_count) * 100.0
		precision = (len(juncs_inport) / len(junc_set)) * 100.0
		print(len(juncs_inport), "/", ref_junc_count, "(" + "{0:.2f}".format(recall) + ") of valid junctions.")
		print(len(juncs_inport), "/", len(junc_set), "(" + "{0:.2f}".format(precision) + ") junctions supported by input.")
		print()


		print("Identifying invalid multi-exonic transcripts ... ", end="")
		invalid_transcripts = collections.defaultdict(list)
		for t, introns in intron_chains.items():
			for index, j in enumerate(introns):
				if not j.key in ref_juncs:
					invalid_transcripts[t].append(str(j.start+1) + "_" + str(j.end+1))
		print(" done.  Found", len(invalid_transcripts), "/", nb_multiexonic, "invalid multi-exonic transcripts.")
		print()


		print("Writing output to", args.output," ... ", end="")
		o = open(args.output, mode='w')

		# Create dict of transcripts with all exons in memory
		with open(args.input[0]) as f:
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

	parser.add_argument("-j", "--junctions", help='''The file containing junctions that should be found in the GTF.  Either this or '-t' must be used.''')
	parser.add_argument("-t", "--transcripts", help='''The file containing transcripts that should be found in the GTF.  Either this or '-j' must be used.''')

	parser.add_argument("-o", "--output", default="junctools.out.gtf",
						help="The filtered or markedup GTF output.  By default we print to stdout.")

	parser.add_argument("mode", help='''GTF operation to apply.  See above for details.  Available options:
 - filter
 - markup
 - compare
''')

	parser.add_argument("input", nargs="+", help="The input GTF file(s) to process")
