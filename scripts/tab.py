#!/usr/bin/env python3

import bed12


class TabEntry:
	id = ""
	refid = ""
	refname = ""
	reflen = 0
	start = 0
	end = 0
	left = 0
	right = 0
	ss1 = ""
	ss2 = ""
	read_strand = "."
	ss_strand = "."
	consensus_strand = "."

	metrics = [29]

	mql = ""
	suspect = False
	pfp = False

	jo = [20]

	def __init__(self):
		self.data = []

	def __str__(self):
		id_parts = [self.id, self.refid, self.refname, self.reflen, self.start, self.end, self.left, self.right,
				self.ss1, self.ss2,
				self.read_strand, self.ss_strand, self.consensus_strand]
		jo_parts = [self.mql, self.suspect, self.pfp]

		chunks = []
		chunks.append("\t".join([str(_) for _ in id_parts]))
		chunks.append("\t".join([str(_) for _ in self.metrics]))
		chunks.append("\t".join([str(_) for _ in jo_parts]))
		chunks.append("\t".join([str(_) for _ in self.jo]))

		return "\t".join(chunks)

	def __key__(self):
		return (self.chrom.encode(), self.start, self.end)

	def __hash__(self):
		return hash(self.__key__())

	def getRaw(self):
		return self.metrics[1]

	def getReliable(self):
		return self.metrics[3]

	def getEntropy(self):
		return self.metrics[10]

	def getEntropyAsStr(self):
		return "{0:.2f}".format(self.getEntropy())

	def getMaxMMES(self):
		return self.metrics[11]

	def getMinHamming(self):
		return min(self.metrics[12], self.metrics[13])

	@property
	def key(self):
		return (self.chrom, self.start, self.end, self.consensus_strand)

	def toExonGFF(self, source="portcullis"):

		entries = []
		parts = [self.chrom, source, "match", self.left + 1, self.right + 1, 0.0, self.strand, ".",
				"ID=" + self.id + ";" +
				"Name=" + self.id + ";" +
				"Note=cov:" + str(self.getRaw()) + "|rel:" + str(self.getReliable()) + "|ent:" + self.getEntropyAsStr() + "|maxmmes:" + str(self.getMaxMMES()) + "|ham:" + str(
					self.getMinHamming()) + ";" +
				"mult=" + str(self.getRaw()) + ";" +
				"grp=" + str(self.id) + ";" +
				"src=E"
				]
		entries.append("\t".join([str(_) for _ in parts]))

		parts = [self.chrom, source, "match_part", self.left + 1, self.start, 0.0, self.strand, ".",
				 "ID=" + self.id + "_left;" +
				 "Parent=" + self.id]
		entries.append("\t".join([str(_) for _ in parts]))

		parts = [self.chrom, source, "match_part", self.end + 2, self.right + 1, 0.0, self.strand, ".",
				 "ID=" + self.id + "_right;" +
				 "Parent=" + self.id]
		entries.append("\t".join([str(_) for _ in parts]))

		return "\n".join(entries)

	def toIntronGFF(self, source="portcullis"):
		parts = [self.chrom, source, "intron", self.start + 1, self.end + 1, self.getRaw(), self.strand, ".",
				#"ID=" + self.id + ";" +
				#"Name=" + self.id + ";" +
				# Removing this as it causes issues with PASA
				#"Note=cov:" + str(self.getRaw()) + "|rel:" + str(self.getReliable()) + "|ent:" + self.getEntropyAsStr() + "|maxmmes:" + str(self.getMaxMMES()) + "|ham:" + str(
				#	self.getMinHamming()) + ";" +
				"mult=" + str(self.getRaw()) + ";" +
				"grp=" + self.id + ";" +
				"src=E"]
		return "\t".join([str(_) for _ in parts])

	@staticmethod
	def metric_names():
		return ["M1-canonical_ss",
				"M2-nb_reads",
				"M3-nb_dist_aln",
				"M4-nb_rel_aln",
				"M5-intron_size",
				"M6-left_anc_size",
				"M7-right_anc_size",
				"M8-max_min_anc",
				"M9-dif_anc",
				"M10-dist_anc",
				"M11-entropy",
				"M12-maxmmes",
				"M13-hamming5p",
				"M14-hamming3p",
				"M15-coverage",
				"M16-uniq_junc",
				"M17-primary_junc",
				"M18-mm_score",
				"M19-mean_mismatches",
				"M20-nb_usrs",
				"M21-nb_msrs",
				"M22-rel2raw",
				"M23-nb_up_juncs",
				"M24-nb_down_juncs",
				"M25-up_aln",
				"M26-down_aln",
				"M27-dist_2_up_junc",
				"M28-dist_2_down_junc",
				"M29-dist_nearest_junc"]

	@staticmethod
	def jo_names():
		return ["JO01",
				"JO02",
				"JO03",
				"JO04",
				"JO05",
				"JO06",
				"JO07",
				"JO08",
				"JO09",
				"JO10",
				"JO11",
				"JO12",
				"JO13",
				"JO14",
				"JO15",
				"JO16",
				"JO17",
				"JO18",
				"JO19",
				"JO20"]

	@staticmethod
	def strand_names():
		return ["read-strand", "ss-strand", "consensus-strand"]

	@staticmethod
	def featureAt(index):
		return TabEntry.features()[index]

	@staticmethod
	def sortedFeatures(indicies):
		f = TabEntry.features()
		s = []
		for i in indicies:
			s.append(f[i])
		return s

	@staticmethod
	def nbMetrics():
		return len(TabEntry.features())

	@staticmethod
	def create_from_tabline(line):
		b = TabEntry()

		parts = line.split("\t")

		assert len(parts >= 56)

		b.id = str(parts[0])
		b.refid = int(parts[1])
		b.refname = parts[2]
		b.reflen = int(parts[3])
		b.start = int(parts[4])
		b.end = int(parts[5])
		b.left = int(parts[6])
		b.right = int(parts[7])
		b.ss1 = parts[8]
		b.ss2 = parts[9]
		b.read_strand = parts[10]
		b.ss_strand = parts[11]
		b.consensus_strand = parts[12]

		b.metrics = parts[13:32]

		b.mql = parts[33]
		b.suspect = bool(parts[34])
		b.pfp = bool(parts[35])

		b.jo = parts[36:56]

		if len(parts > 56):
			i=0

		return b


def loadtab(tabfile):
	bed = list()
	tab = list()
	with open(tabfile) as f:
		# Skip header
		f.readline()

		for line in f:
			line.strip()
			if len(line) > 1:
				bed.append(bed12.BedEntry.create_from_tabline(line, False, False))
				tab.append(TabEntry.create_from_tabline(line))

	return bed, tab


def loadtabasset(tabfile):
	tab = set()
	with open(tabfile) as f:
		# Skip header
		f.readline()

		for line in f:
			line.strip()
			if len(line) > 1:
				t = TabEntry.create_from_tabline(line)
				tab.add(t.key)

	return tab


def filtertab(filepath, outfile, tab_set, mode):
	o = open(outfile, 'w')

	index = 0;
	with open(filepath) as f:

		o.write(f.readline())

		for line in f:

			line = line.strip()
			if line != "":
				key = TabEntry.create_from_tabline(line).key
				if mode == 0:
					if key in tab_set:
						o.write(line + "\n")
				elif mode == 1:
					if key not in tab_set:
						o.write(line + "\n")

	o.close()
