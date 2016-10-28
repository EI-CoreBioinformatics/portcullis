#!/usr/bin/env python3

import os
import abc

__author__ = 'maplesod'


class Junction(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, use_strand=True):
		self._use_strand = use_strand
		self.refseq = ""
		self.start = 0
		self.end = 0
		self.strand = "."
		self.score = 0.0
		self.id = ""

	def __str__(self):
		line = [self.id, self.refseq, self.start, self.end, self.strand, self.score]
		return "\t".join([str(_) for _ in line])

	def __cmp__(self, other):
		if hasattr(other, 'refseq') and hasattr(other, 'start') and hasattr(other, 'end'):
			if self.__lt__(other):
				return 1
			elif self.__gt__(other):
				return -1
			else:
				return 0

	def __key__(self):
		return (self.refseq.encode(), self.start, self.end, self.strand if self._use_strand else None)

	def __hash__(self):
		return hash(self.__key__())

	@property
	def key(self):
		return self.__key__()

	def __lt__(self, other):
		if self.refseq.__lt__(other.refseq):
			return True
		else:
			if self.refseq.__gt__(other.refseq):
				return False
			else:
				# The same chrom
				if self.start < other.start:
					return True
				else:
					if self.start > other.start:
						return False
					else:
						# Also the same start
						if self.end < other.end:
							return True
						else:
							if self.end > other.end:
								return False
							else:
								# Same!
								return False

	def __eq__(self, other):
		return not self < other and not other < self

	def __ne__(self, other):
		return self < other or other < self

	def __gt__(self, other):
		return other < self

	def __ge__(self, other):
		return not self < other

	def __le__(self, other):
		return not other < self

	@abc.abstractmethod
	def parse_line(line, fullparse=True):
		pass

	@abc.abstractmethod
	def accepts_ext(ext):
		pass

	@abc.abstractmethod
	def file_header(self, description=None):
		pass

	@staticmethod
	def createDict(self, filepath):

		index = 0
		items = {}

		filename, ext = os.path.splitext(filepath)

		with open(filepath) as f:
			for line in f:
				junc = create_junction(ext)
				res = junc.parse_line(line, fullparse=False)
				if res == None:
					continue
				items[junc.key] = line
				index += 1
		if len(items) != index:
			print("duplicated items in file " + filepath)
		return items

	@staticmethod
	def createSet(self, filepath, keyonly=False):

		index = 0
		items = set()

		filename, ext = os.path.splitext(filepath)

		with open(filepath) as f:
			for line in f:
				junc = create_junction(ext)
				res = junc.parse_line(line, fullparse=not keyonly)
				if res == None:
					continue
				items.add(junc.key)
				index += 1
		if len(items) != index:
			print("duplicated items in file " + filepath)
		return items

	@staticmethod
	def sort(junctions):
		# Sort by
		junctions.sort(key=lambda x: x.strand)
		junctions.sort(key=lambda x: x.end)
		junctions.sort(key=lambda x: x.start)
		junctions.sort(key=lambda x: x.refseq)

	@staticmethod
	def reindex(junctions, prefix="", start=0):
		index=start
		for j in junctions:
			j.id = prefix + str(index)
			index += 1



class ExonJunction(Junction):
	__metaclass__ = abc.ABCMeta

	def __init__(self, use_strand=True):
		Junction.__init__(self, use_strand=use_strand)
		self.left = 0
		self.right = 0

	def __str__(self):
		line = [super.__str__(), self.left, self.right]
		return "\t".join([str(_) for _ in line])


class Bed6Junction(Junction):
	def __init__(self, use_strand=True):
		Junction.__init__(self, use_strand=use_strand)

	def __str__(self):
		line = [self.refseq, self.start, self.end, self.id, self.score,
				self.strand if self.strand else "."]
		return "\t".join([str(_) for _ in line])

	def accepts_ext(ext):
		return ext == ".bed"

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		# Handle header or blank lines
		if (len(parts) != 6):
			return None

		self.refseq = parts[0]
		self.strand = parts[5]
		self.start = int(parts[1])
		self.end = int(parts[2])

		if fullparse:
			self.id = parts[3]
			self.score = float(parts[4])

		return self

	def file_header(self, description=None):
		d = ""
		if not description == None and not description == "":
			d = "description=\"{0}\"".format(description)

		return "track name=\"junctions\"" + d


class Bed12Junction(ExonJunction):
	def __init__(self, use_strand=True, tophat=False):
		ExonJunction.__init__(self, use_strand=use_strand)
		self.tophat = tophat
		self.name = ""
		self.red = 0
		self.green = 0
		self.blue = 0
		self.block_count = 0
		self.block_sizes = []
		self.block_starts = []

	def __str__(self):

		rgb = ",".join([str(_) for _ in (self.red, self.green, self.blue)])
		assert len(self.block_sizes) == len(self.block_starts) == self.block_count, (self.block_count,
																					 len(self.block_sizes),
																					 len(self.block_starts))
		bsizes = ",".join([str(_) for _ in self.block_sizes])
		bstarts = ",".join([str(_) for _ in self.block_starts])

		line = [self.refseq, self.left, self.right, self.name, self.score,
				self.strand if self.strand else ".",
				self.start, self.end,
				rgb,
				self.block_count,
				bsizes,
				bstarts
				]
		return "\t".join([str(_) for _ in line])

	def toIntronStyle(self, bed6=False):

		rgb = ",".join([str(_) for _ in (self.red, self.green, self.blue)])
		assert len(self.block_sizes) == len(self.block_starts) == self.block_count, (self.block_count,
																					 len(self.block_sizes),
																					 len(self.block_starts))
		line = ""
		if not bed6:
			line = [self.refseq, self.start, self.end, self.name, self.score,
					self.strand if self.strand else ".",
					self.start, self.end,
					rgb,
					2,
					"0,0",
					"0,0"
					]
		else:
			line = [self.refseq, self.start, self.end, self.name, self.score,
					self.strand if self.strand else "."]
		return "\t".join([str(_) for _ in line])

	def accepts_ext(ext):
		return ext == ".bed"

	def file_header(self, description=None):
		d = ""
		if not description == None and not description == "":
			d = "description=\"{0}\"".format(description)

		return "track name=\"junctions\"" + d

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		# Handle header or blank lines
		if (len(parts) != 12):
			return None

		self.refseq = parts[0]
		self.strand = parts[5]
		self.start = int(parts[6])
		self.end = int(parts[7])

		if fullparse:
			self.left = int(parts[1])
			self.right = int(parts[2])
			self.name = parts[3]
			self.score = float(parts[4])

			c_parts = parts[8].split(",")
			self.red = int(c_parts[0])
			self.green = int(c_parts[1])
			self.blue = int(c_parts[2])
			self.block_count = int(parts[9])

			self.block_sizes = [int(_) for _ in parts[10].split(",")]
			self.block_starts = [int(_) for _ in parts[11].rstrip().split(",")]

			if self.tophat:
				self.start += self.block_sizes[0]
				self.end -= self.block_sizes[1]

			assert len(self.block_sizes) == len(self.block_starts) == self.block_count, (line,
																						 self.block_count,
																						 self.block_sizes,
																						 self.block_starts)
		return self

	@staticmethod
	def create_from_tabline(key, use_strand=True):

		b = Bed12Junction(use_strand=use_strand)

		parts = key.split("\t")

		b.chrom = parts[2]
		b.start = int(parts[6])
		b.end = int(parts[7]) + 1
		b.strand = parts[12]
		b.score = float(parts[14])
		b.thick_start = int(parts[4])
		b.thick_end = int(parts[5]) + 1
		b.name = "junc"

		b.red = 255
		b.green = 0
		b.blue = 0
		b.block_count = 2

		b.block_sizes = [(b.thick_start - b.start + 1), (b.end - b.thick_end + 1)]
		# for bp in bsize_parts:
		#     b.block_sizes.append(int(bp))

		b.block_starts = [0, b.thick_end - b.start]

		# for bp in bstart_parts:
		#     b.block_starts.append(int(bp))
		assert len(b.block_sizes) == len(b.block_starts) == b.block_count, (key,
																			b.block_count,
																			b.block_sizes,
																			b.block_starts)

		return b


class TabJunction(ExonJunction):
	def __init__(self, use_strand=True):
		ExonJunction.__init__(self, use_strand=use_strand)
		self.id = ""
		self.refid = ""
		self.reflen = 0
		self.ss1 = ""
		self.ss2 = ""
		self.read_strand = "."
		self.ss_strand = "."

		self.metrics = [len(TabJunction.metric_names())]

		self.mql = ""
		self.suspect = False
		self.pfp = False

		self.jo = [len(TabJunction.jo_names())]

	def __str__(self):
		id_parts = [self.id, self.refid, self.refseq, self.reflen, self.start, self.end, self.left, self.right,
					self.ss1, self.ss2,
					self.read_strand, self.ss_strand, self.strand]
		jo_parts = [self.mql, self.suspect, self.pfp]

		chunks = []
		chunks.append("\t".join([str(_) for _ in id_parts]))
		chunks.append("\t".join([str(_) for _ in self.metrics]))
		chunks.append("\t".join([str(_) for _ in jo_parts]))
		chunks.append("\t".join([str(_) for _ in self.jo]))

		return "\t".join(chunks)

	def accepts_ext(ext):
		return ext == ".tab"

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

	def toExonGFF(self, source="portcullis"):

		entries = []
		parts = [self.chrom, source, "match", self.left + 1, self.right + 1, 0.0, self.strand, ".",
				 "ID=" + self.id + ";" +
				 "Name=" + self.id + ";" +
				 "Note=cov:" + str(self.getRaw()) + "|rel:" + str(
					 self.getReliable()) + "|ent:" + self.getEntropyAsStr() + "|maxmmes:" + str(
					 self.getMaxMMES()) + "|ham:" + str(
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
				 # "ID=" + self.id + ";" +
				 # "Name=" + self.id + ";" +
				 # Removing this as it causes issues with PASA
				 # "Note=cov:" + str(self.getRaw()) + "|rel:" + str(self.getReliable()) + "|ent:" + self.getEntropyAsStr() + "|maxmmes:" + str(self.getMaxMMES()) + "|ham:" + str(
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

	def file_header(self, description=""):
		chunks = []
		chunks.append("\t".join(["index", "refid", "refname", "reflen", "start", "end", "left", "right", "ss1", "ss2"]))
		chunks.append("\t".join(TabJunction.strand_names()))
		chunks.append("\t".join(TabJunction.metric_names()))
		chunks.append("\t".join(["MQL", "Suspect", "PFP"]))
		chunks.append("\t".join(TabJunction.jo_names()))

		return "\t".join(chunks)

	@staticmethod
	def metric_at(index):
		return TabJunction.mefeatures()[index]

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if len(parts) < 56:
			return None
		elif parts[0] == "index":
			return None

		self.refseq = parts[2]
		self.start = int(parts[4])
		self.end = int(parts[5])
		self.strand = parts[12]

		if fullparse:
			self.id = str(parts[0])
			self.refid = int(parts[1])
			self.reflen = int(parts[3])
			self.left = int(parts[6])
			self.right = int(parts[7])
			self.ss1 = parts[8]
			self.ss2 = parts[9]
			self.read_strand = parts[10]
			self.ss_strand = parts[11]

			self.metrics = parts[13:32]

			self.mql = parts[33]
			self.suspect = bool(parts[34])
			self.pfp = bool(parts[35])

			self.jo = parts[36:56]

			if len(parts) > 56:
				i = 0

		return self


def saveList(filepath, list, description):
	o = open(filepath, 'w')

	o.write("track name=\"junctions\" description=\"" + description + "\"\n")

	for b in list:
		print(b, file=o)
	o.close()


def filterbed(filepath, refset, mode, usestrand, tophat, outfile):
	o = open(outfile, 'w')

	o.write("track name=\"junctions\"\n")

	index = 0
	with open(filepath) as i:

		i.readline()  # Skip header
		for line in i:

			key = Bed12Junction(use_strand=usestrand, tophat=tophat).parse_line(line, fullparse=False).key

			if mode == 0:
				if key in refset:
					o.write(line)
			elif mode == 1:
				if key not in refset:
					o.write(line)

			index += 1
	o.close()


def filtertab(filepath, outfile, bed, mode, usestrand):
	o = open(outfile, 'w')

	index = 0;
	with open(filepath) as f:

		o.write(f.readline())

		for line in f:

			line = line.strip()
			if line != "":
				key = TabJunction(use_strand=usestrand).parse_line(line, fullparse=False).key
				if mode == 0:
					if key in bed:
						o.write(line + "\n")
				elif mode == 1:
					if key not in bed:
						o.write(line + "\n")

	o.close()


def create_junction(ext, use_strand=True):
	for cls in Junction.__subclasses__():
		if cls.accepts_ext(ext):
			return cls(use_strand)
	raise ValueError


def create_exon_junction(ext, use_strand=True):
	for cls in ExonJunction.__subclasses__():
		if cls.accepts_ext(ext):
			return cls(use_strand)
	raise ValueError
