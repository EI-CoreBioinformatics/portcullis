import os
import argparse
import bed12


def set2list(junc_set, name_prefix):
	new_bed = list()
	for e in junc_set:
		new_bed.append(bed12.BedEntry.create_from_key(e))
	new_bed.sort()
	i = 1
	for b in new_bed:
		b.name = name_prefix + "_junc_" + str(i)
		i += 1

	return new_bed


def main():
	parser = argparse.ArgumentParser("Script to create the Venn Plots from BED files")
	parser.add_argument("input", help="The directory containing BED files from pipeline")
	parser.add_argument("-o", "--output", required=True, help="The output bed files prefix")
	parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
	args = parser.parse_args()

	ref_bed = bed12.loadbed(args.reference, False, False)
	print("Loaded Reference BED file.  # junctions: " + str(len(ref_bed)))

	# Load all bed files
	bed_data = {}
	aligners = set()
	reads = set()
	junc_analysers = set()
	for bed_file in os.listdir(args.input):
		if bed_file.endswith("-real-all.bed"):
			bed_base = os.path.splitext(bed_file)[0]
			bed_data[bed_base] = bed12.loadbed(args.input + "/" + bed_file, False, False)
			parts = bed_base.split('-')
			aligners.add(parts[0])
			reads.add(parts[1])
			junc_analysers.add(parts[2])
			print("Loaded: " + bed_file + "; # junctions: " + str(len(bed_data[bed_base])))

	print("Found these aligners: " + ', '.join(aligners))
	print("Found these reads: " + ', '.join(reads))
	print("Found these junction analysis tools: " + ', '.join(junc_analysers))

	# Filtering reference
	new_ref = set()
	extra = set()
	for a in aligners:
		new_ref = new_ref.union(ref_bed.intersection(bed_data[a + "-real-all"]))
		extra = extra.union(bed_data[a + "-real-all"] - ref_bed)

	print("New reference contains " + str(len(new_ref)) + " junctions from original reference")

	sv_bed = set2list(new_ref, "simvar")

	# Output new bed file to disk
	with open(args.output + ".sim_var.bed", "w") as bed_sv_out:
		print("track name=\"junctions\"", file=bed_sv_out)
		for b in sv_bed:
			print(b, file=bed_sv_out)
	print("Saved: " + args.output + ".sim_var.bed")
	print()
	print("Found " + str(len(extra)) + " potential junctions outside reference")

	extra_special = set()
	for e in extra:
		found = 0
		for a in aligners:
			if e in bed_data[a + "-real-all"]:
				found += 1

		if found >= 2:
			extra_special.add(e)

	print("Found " + str(len(extra_special)) + " junctions outside reference but present in at least 2 alignments")

	new_ref = new_ref.union(extra_special)

	print("New reference contains " + str(len(new_ref)) + " junctions")

	new_bed = set2list(new_ref, "real")

	# Output new bed file to disk
	with open(args.output + ".real.bed", "w") as bed_out:
		print("track name=\"junctions\"", file=bed_out)
		for b in new_bed:
			print(b, file=bed_out)

	print("Saved: " + args.output + ".real.bed")


main()
