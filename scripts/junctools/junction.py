#!/usr/bin/env python3

import abc
import copy
import os
from enum import Enum, unique

__author__ = 'maplesod'


@unique
class JuncFactory(Enum):
	PORTCULLIS = 1
	BED = 2
	EBED = 3
	TBED = 4
	IBED = 5
	BED6 = 6
	GFF = 7
	EGFF = 8
	IGFF = 9
	GTF = 10
	STAR = 11
	HISAT = 12
	FINESPLICE = 13
	SOAPSPLICE = 14
	SPANKI = 15
	TRUESIGHT = 16
	MAPSPLICE = 17

	def isStreamable(self):
		return False if self == JuncFactory.GFF or self == JuncFactory.EGFF or self == JuncFactory.GTF else True

	def exon_based(self):
		return True if self.PORTCULLIS or self.EBED or self.TBED or self.EGFF else False

	def isBed(self):
		return self.value >= JuncFactory.BED.value and self.value <= JuncFactory.BED6.value

	def isGFF(self):
		return self.value >= JuncFactory.GFF.value and self.value <= JuncFactory.IGFF.value

	@staticmethod
	def create_from_file(filepath, use_strand=True, junc_to_copy=None):
		'''
		Creates an instance of a junction that
		:param filepath:
		:param use_strand:
		:return:
		'''
		filename, ext = os.path.splitext(filepath)
		return JuncFactory.create_from_ext(ext, use_strand=use_strand, junc_to_copy=junc_to_copy)

	@staticmethod
	def create_from_ext(ext, use_strand=True, junc_to_copy=None):
		for cls in ExonJunction.__subclasses__():
			if cls.accepts_ext(ext):
				return cls(use_strand=use_strand, junc_to_copy=junc_to_copy)
		raise ValueError("No junction can be created for " + ext + " files")

	@staticmethod
	def create_from_enum(type, use_strand=True, junc_to_copy=None):
		if type == JuncFactory.PORTCULLIS:
			return TabJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type.isBed():
			return BedJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type.isGFF():
			return GFFJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type == JuncFactory.STAR:
			return StarJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type == JuncFactory.HISAT:
			return HisatJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type == JuncFactory.FINESPLICE:
			return FinespliceJunction(junc_to_copy=junc_to_copy)
		elif type == JuncFactory.SOAPSPLICE:
			return SoapspliceJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type == JuncFactory.SPANKI:
			return SpankiJunction(use_strand=use_strand, junc_to_copy=junc_to_copy)
		elif type == JuncFactory.TRUESIGHT:
			return TruesightJunction(junc_to_copy=junc_to_copy)
		elif type == JuncFactory.MAPSPLICE:
			return MapspliceJunction(junc_to_copy=junc_to_copy)
		raise "Unknown type"


class Junction(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, use_strand=True, junc_to_copy=None):

		if junc_to_copy:
			self._use_strand = junc_to_copy._use_strand
			self.refseq = junc_to_copy.refseq
			self.start = junc_to_copy.start
			self.end = junc_to_copy.end
			self.strand = junc_to_copy.strand
			self.score = junc_to_copy.score
			self.id = junc_to_copy.id
			self.canonical = junc_to_copy.canonical
		else:
			self._use_strand = use_strand
			self.refseq = ""
			self.start = 0
			self.end = 0
			self.strand = "."
			self.score = 0.0
			self.id = ""
			self.canonical = ""

	def __str__(self):
		line = [self.id, self.refseq, self.start, self.end, self.strand, self.score, self.canonical]
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

	def startSplicesiteKey(self):
		return (self.refseq.encode(), self.start, self.strand if self._use_strand else None)

	def endSplicesiteKey(self):
		return (self.refseq.encode(), self.end, self.strand if self._use_strand else None)

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

	def size(self):
		return self.end - self.start + 1

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
	def createDict(filepath, use_strand=True, fullparse=False):

		index = 0
		items = {}

		with open(filepath) as f:
			for line in f:
				junc = JuncFactory.create_from_file(filepath, use_strand=use_strand).parse_line(line,
																								fullparse=fullparse)
				if junc:
					items[junc.key] = line if not fullparse else junc
					index += 1
		return items, index

	@staticmethod
	def createJuncSet(filepath, use_strand=True, fullparse=False):

		index = 0
		items = set()

		with open(filepath) as f:
			for line in f:
				junc = JuncFactory.create_from_file(filepath, use_strand=use_strand).parse_line(line,
																								fullparse=fullparse)
				if junc:
					items.add(junc.key)
					index += 1
		return items, index

	@staticmethod
	def createMarkedupJuncSets(filepath, labelpath, use_strand=True, fullparse=False):
		# Load reference and labels, divide into tp and tn
		rp = set()
		rn = set()

		labs = open(labelpath, "r")
		juncs = open(filepath, "r")

		# Assume we have a header in the juncs file
		header = juncs.readline()

		# Now we should have the same number of lines in both files
		line = 1
		while 1:
			ref_line = juncs.readline()
			lab_line = labs.readline()
			if not ref_line and not lab_line: break
			if (lab_line and not ref_line) or (not lab_line and ref_line): print(
				"ERROR: reference file and labels file have a different number of entries.", file=sys.stderr); exit(1)
			ref_line = ref_line.strip()
			lab_line = lab_line.strip()
			junc = JuncFactory.create_from_file(filepath, use_strand=use_strand).parse_line(line, fullparse=fullparse)
			if lab_line == "1":
				rp.add(junc)
			elif lab_line == "0":
				rn.add(junc)
			else:
				raise ValueError(
					"ERROR: Label file contains an entry that is not either \"0\" or \"1\" at line:" + line)
			line += 1
		labs.close()
		juncs.close()

		return rp, rn, line - 1

	@staticmethod
	def createSpliceSiteSet(filepath, use_strand=True, fullparse=False):

		items = {}

		with open(filepath) as f:
			for line in f:
				junc = JuncFactory.create_from_file(filepath, use_strand=use_strand).parse_line(line,
																								fullparse=fullparse)

				if junc:

					key1 = junc.startSplicesiteKey()
					key2 = junc.endSplicesiteKey()

					if key1 in items:
						items[key1] += 1
					else:
						items[key1] = 1

					if key2 in items:
						items[key2] += 1
					else:
						items[key2] = 1
		return items

	@staticmethod
	def save(filepath, junctions, description, bed6=False):
		o = open(filepath, 'w')

		junc = Junction.create(filepath)
		o.write(junc.file_header(description=description))

		if bed6 and not junc.accepts_ext(".bed"):
			raise "Requested BED6 output but output file does not have a bed extension"

		for b in junctions:
			if bed6:
				o.write(b.toIntronStyle(bed6=True))
			else:
				o.write(b)
		o.close()

	@staticmethod
	def sort(junctions):
		# Sort by
		junctions.sort(key=lambda x: x.strand)
		junctions.sort(key=lambda x: x.end)
		junctions.sort(key=lambda x: x.start)
		junctions.sort(key=lambda x: x.refseq)

	@staticmethod
	def reindex(junctions, prefix="", start=0):
		index = start
		for j in junctions:
			j.id = prefix + str(index)
			index += 1


class ExonJunction(Junction):
	__metaclass__ = abc.ABCMeta

	def __init__(self, use_strand=True, junc_to_copy=None):
		Junction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

		if junc_to_copy and issubclass(type(junc_to_copy), ExonJunction):
			self.left = junc_to_copy.left
			self.right = junc_to_copy.right
		else:
			self.left = 0
			self.right = 0

	def __str__(self):
		line = [super.__str__(), self.left, self.right]
		return "\t".join([str(_) for _ in line])


class BedJunction(ExonJunction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		ExonJunction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

		if junc_to_copy and type(junc_to_copy) is BedJunction:
			self.red = junc_to_copy.red
			self.green = junc_to_copy.green
			self.blue = junc_to_copy.blue
			self.style = junc_to_copy.style
		else:
			self.red = 255
			self.green = 0
			self.blue = 0
			self.style = JuncFactory.IBED

	def __str__(self):

		line = ""

		scr = "{0:.3f}".format(self.score)

		if self.style == JuncFactory.BED6:

			line = [self.refseq, self.start, self.end + 1, self.id, scr,
					self.strand if self.strand else "."]

		else:

			rgb = ",".join([str(_) for _ in (self.red, self.green, self.blue)])

			if self.style == JuncFactory.IBED:
				line = [self.refseq, self.start, self.end + 1, self.id, scr,
						self.strand if self.strand else ".",
						self.start, self.end + 1,
						rgb,
						2,
						"0,0",
						"0,0"
						]
			else:

				block_count = 2
				bsizes = ",".join([self.start - self.left, self.right - self.end])
				bstarts = ",".join([0, self.end - self.left + 1])


				if self.style == JuncFactory.EBED:

					line = [self.refseq, self.left, self.right + 1, self.id, scr,
							self.strand if self.strand else ".",
							self.start, self.end + 1,
							rgb,
							block_count,
							bsizes,
							bstarts
							]
				elif self.style == JuncFactory.TBED:

					line = [self.refseq, self.left, self.right + 1, self.id, scr,
							self.strand if self.strand else ".",
							self.left, self.right + 1,
							rgb,
							block_count,
							bsizes,
							bstarts
							]

		return "\t".join([str(_) for _ in line])

	def accepts_ext(ext):
		return ext == ".bed" or ext == ".ibed" or ext == ".ebed" or ext == ".bed6"

	def file_header(self, description=None):
		d = ""
		if not description == None and not description == "":
			d = "description=\"{0}\"".format(description)

		return "track name=\"junctions\"" + d

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		# Handle header or blank lines
		if (len(parts) != 6 and len(parts) != 12):
			return None

		self.style = JuncFactory.BED6 if len(parts) == 6 else JuncFactory.IBED

		self.refseq = parts[0]
		self.strand = parts[5]
		self.start = int(parts[1]) if self.style == JuncFactory.BED6 else int(parts[6])
		self.end = int(parts[2]) - 1 if self.style == JuncFactory.BED6 else int(parts[7]) - 1

		if fullparse:
			self.id = parts[3]
			self.score = float(parts[4])

			if self.style != JuncFactory.BED6:
				self.left = int(parts[1])
				self.right = int(parts[2]) - 1

				c_parts = parts[8].split(",")
				self.red = int(c_parts[0])
				self.green = int(c_parts[1])
				self.blue = int(c_parts[2])
				self.block_count = int(parts[9])

				block_sizes = [int(_) for _ in parts[10].split(",")]

				# Check if this looks like a tophat style junction and if so bring it into out style
				if self.start == self.left and block_sizes[0] != 0:
					self.style = JuncFactory.TBED
					self.start += block_sizes[0]
					self.end -= block_sizes[1]
				elif self.start != self.left:
					self.style = JuncFactory.EBED

		return self


class GFFJunction(ExonJunction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		ExonJunction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

		self.style = JuncFactory.IGFF
		self.source = "portcullis"
		self.feature = "intron"
		self.frame = "."
		self.attrs = []
		self.note = ""
		self.raw = 0

		if junc_to_copy:
			if type(junc_to_copy) is TabJunction:
				self.note = "Note=can:" + junc_to_copy.getSSType() + "|cov:" + str(junc_to_copy.getRaw()) + "|rel:" + str(
					junc_to_copy.getReliable()) + "|ent:" + junc_to_copy.getEntropyAsStr() + "|maxmmes:" + str(
					junc_to_copy.getMaxMMES()) + "|ham:" + str(
					junc_to_copy.getMinHamming()) + ";"
				self.id = junc_to_copy.id
				self.score = junc_to_copy.getScore()
				self.raw = junc_to_copy.getRaw()
			elif type(junc_to_copy) is GFFJunction:
				self.note = junc_to_copy.note
				self.score = junc_to_copy.score
				self.raw = junc_to_copy.raw
				self.source = junc_to_copy.source
				self.feature = junc_to_copy.feature
				self.frame = junc_to_copy.frame
				self.style = junc_to_copy.style

	def __str__(self):

		if self.style == JuncFactory.EGFF:
			entries = []
			parts = [self.refseq, self.source, "match", self.left + 1, self.right + 1, self.score, self.strand,
					 self.frame,
					 "ID=" + self.id + ";" +
					 "Name=" + self.id + ";" +
					 self.note +
					 "mult=" + str(self.raw) + ";" +
					 "grp=" + str(self.id) + ";" +
					 "src=E"
					 ]
			entries.append("\t".join([str(_) for _ in parts]))

			parts = [self.refseq, self.source, "match_part", self.left + 1, self.start, 0.0, self.strand, self.frame,
					 "ID=" + self.id + "_left;" +
					 "Parent=" + self.id]
			entries.append("\t".join([str(_) for _ in parts]))

			parts = [self.refseq, self.source, "match_part", self.end + 2, self.right + 1, 0.0, self.strand, self.frame,
					 "ID=" + self.id + "_right;" +
					 "Parent=" + self.id]
			entries.append("\t".join([str(_) for _ in parts]))

			return "\n".join(entries)
		else:
			parts = [self.refseq, self.source, self.feature, self.start + 1, self.end + 1, self.score, self.strand,
					 self.frame,
					 # "ID=" + self.id + ";" +
					 # "Name=" + self.id + ";" +
					 # Removing this as it causes issues with PASA
					 # "Note=cov:" + str(self.getRaw()) + "|rel:" + str(self.getReliable()) + "|ent:" + self.getEntropyAsStr() + "|maxmmes:" + str(self.getMaxMMES()) + "|ham:" + str(
					 #	self.getMinHamming()) + ";" +
					 "mult=" + str(self.raw) + ";" +
					 "grp=" + self.id + ";" +
					 "src=E"]
		return "\t".join([str(_) for _ in parts])

	def parse_line(self, line, fullparse=True):

		if line.startswith("#"):
			return None

		parts = line.split("\t")

		if len(parts) <= 1:
			return None

		if len(parts) != 9 and len(parts) > 1:
			msg = "Unexpected number of columns in GFF file.  Expected 9, found " + str(len(parts))
			raise ValueError(msg)

		if parts[2] != "intron":
			return None

		self.refseq = parts[0]
		self.start = int(parts[3]) - 1
		self.end = int(parts[4]) - 1
		self.strand = parts[6]

		if fullparse:
			self.source = parts[1]
			self.feature = parts[2]
			self.score = float(parts[5])
			self.frame = parts[7]
			self.attrs = parts[8].split(";")
			for a in self.attrs:
				p = a.split("=")
				nam = p[0]
				val = p[1]
				if nam.startswith("ID"):
					self.id = val
				elif nam.startswith("mult"):
					self.raw = int(val)
				elif nam.startswith("Note"):
					self.note = val

		return self


class TabJunction(ExonJunction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		ExonJunction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

		if junc_to_copy and type(junc_to_copy) is TabJunction:
			self.refid = junc_to_copy.refid
			self.reflen = junc_to_copy.reflen
			self.ss1 = junc_to_copy.ss1
			self.ss2 = junc_to_copy.ss2
			self.read_strand = junc_to_copy.read_strand
			self.ss_strand = junc_to_copy.ss_strand

			self.metrics = copy.deepcopy(junc_to_copy.metrics)
			self.jo = copy.deepcopy(junc_to_copy.jo)
		else:
			self.refid = ""
			self.reflen = 0
			self.ss1 = ""
			self.ss2 = ""
			self.read_strand = "."
			self.ss_strand = "."

			self.metrics = [len(TabJunction.metric_names())]
			self.jo = [len(TabJunction.jo_names())]

	def __str__(self):
		id_parts = [self.id, self.refid, self.refseq, self.reflen, self.start, self.end, self.size(), self.left, self.right,
					self.read_strand, self.ss_strand, self.strand, self.ss1, self.ss2]

		chunks = []
		chunks.append("\t".join([str(_) for _ in id_parts]))
		chunks.append("\t".join([str(_) for _ in self.metrics]))
		chunks.append("\t".join([str(_) for _ in self.jo]))

		return "\t".join(chunks)

	def accepts_ext(ext):
		return ext == ".tab"

	def getSSType(self):
		return self.metrics[0]

	def getRaw(self):
		return int(self.metrics[4])

	def getScore(self):
		return float(self.metrics[1])

	def getScoreAsStr(self):
		return "{0:.2f}".format(self.getScore())

	def getReliable(self):
		return int(self.metrics[12])

	def getEntropy(self):
		return float(self.metrics[14])

	def getEntropyAsStr(self):
		return "{0:.2f}".format(self.getEntropy())

	def getMaxMMES(self):
		return int(self.metrics[18])

	def getMinHamming(self):
		return min(int(self.metrics[20]), int(self.metrics[21]))

	@staticmethod
	def metric_names():
		return ["canonical_ss",
				"score",
				"suspicious",
				"pfp",
				"nb_raw_aln",
				"nb_dist_aln",
				"nb_us_aln",
				"nb_ms_aln",
				"nb_um_aln",
				"nb_mm_aln",
				"nb_bpp_aln",
				"nb_ppp_aln",
				"nb_rel_aln",
				"rel2raw",
				"entropy",
				"mean_mismatches",
				"mean_readlen",
				"max_min_anc",
				"maxmmes",
				"intron_score",
				"hamming5p",
				"hamming3p",
				"coding",
				"pws",
				"splice_sig",
				"uniq_junc",
				"primary_junc",
				"nb_up_juncs",
				"nb_down_juncs",
				"dist_2_up_junc",
				"dist_2_down_junc",
				"dist_nearest_junc",
				"mm_score",
				"coverage",
				"up_aln",
				"down_aln"]

	@staticmethod
	def jo_names():
		return ["JAD01",
				"JAD02",
				"JAD03",
				"JAD04",
				"JAD05",
				"JAD06",
				"JAD07",
				"JAD08",
				"JAD09",
				"JAD10",
				"JAD11",
				"JAD12",
				"JAD13",
				"JAD14",
				"JAD15",
				"JAD16",
				"JAD17",
				"JAD18",
				"JAD19",
				"JAD20"]

	@staticmethod
	def strand_names():
		return ["read-strand", "ss-strand", "consensus-strand"]

	def file_header(self, description=""):
		chunks = []
		chunks.append("\t".join(["index", "refid", "refname", "reflen", "start", "end", "size", "left", "right"]))
		chunks.append("\t".join(TabJunction.strand_names()))
		chunks.append("\t".join(["ss1", "ss2"]))
		chunks.append("\t".join(TabJunction.metric_names()))
		chunks.append("\t".join(TabJunction.jo_names()))

		return "\t".join(chunks)

	@staticmethod
	def metric_at(index):
		return TabJunction.mefeatures()[index]

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if parts[0] == "index" or len(parts) <= 1:
			return None

		if len(parts) != 70 and len(parts) > 1:
			msg = "Unexpected number of columns in TAB file.  Expected 70, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[2]
		self.start = int(parts[4])
		self.end = int(parts[5])
		self.strand = parts[11]

		if fullparse:
			self.id = str(parts[0])
			self.refid = int(parts[1])
			self.reflen = int(parts[3])
			self.left = int(parts[7])
			self.right = int(parts[8])
			self.ss1 = parts[12]
			self.ss2 = parts[13]
			self.read_strand = parts[9]
			self.ss_strand = parts[10]

			jostart = 14+len(TabJunction.metric_names())
			self.metrics = parts[14:jostart]

			endpart = jostart+len(TabJunction.jo_names())
			self.jo = parts[jostart:endpart]

			if len(parts) > 56:
				i = 0

		return self


class StarJunction(Junction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		Junction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

		if junc_to_copy and type(junc_to_copy) is StarJunction:
			self.motif = junc_to_copy.motif
			self.annotated = junc_to_copy.annotated
			self.mm = junc_to_copy.mm
			self.overhang = junc_to_copy.overhang

		else:
			self.motif = ""
			self.annotated = 0
			self.mm = 0
			self.overhang = 0

	def __str__(self):
		st = 1 if self.strand == "+" else 2 if self.strand == "-" else 0
		line = [self.refseq, self.start + 1, self.end + 1, st, self.motif, self.annotated, self.score, self.mm,
				self.overhang]
		return "\t".join([str(_) for _ in line])

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if len(parts) <= 1:
			return None

		if len(parts) != 9 and len(parts) > 1:
			msg = "Unexpected number of columns in STAR junction file.  Expected 9, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[0]
		self.start = int(parts[1]) - 1
		self.end = int(parts[2]) - 1
		self.strand = "+" if parts[3] == "1" else "-" if parts[3] == "2" else "."

		if fullparse:
			self.motif = parts[4]
			self.annotated = int(parts[5])
			self.score = int(parts[6])
			self.mm = int(parts[7])
			self.overhang = int(parts[8])

		return self


class HisatJunction(Junction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		Junction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

	def __str__(self):
		line = [self.refseq, self.start - 1, self.end + 1, self.strand]
		return "\t".join([str(_) for _ in line])

	def parse_line(self, line, fullparse=True):

		parts = line.strip().split("\t")

		if len(parts) <= 1:
			return None

		if len(parts) != 4 and len(parts) > 1:
			msg = "Unexpected number of columns in Hisat junction file.  Expected 4, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[0]
		self.start = int(parts[1]) + 1
		self.end = int(parts[2]) - 1
		self.strand = parts[3]

		return self


class FinespliceJunction(Junction):
	def __init__(self, junc_to_copy=None):
		Junction.__init__(self, use_strand=False, junc_to_copy=junc_to_copy)

		if junc_to_copy and type(junc_to_copy) is FinespliceJunction:
			self.unique = junc_to_copy.unique
			self.rescued = junc_to_copy.rescued
		else:
			self.unique = 0
			self.rescued = 0

	def __str__(self):
		line = [self.refseq, self.start, self.end + 1, self.score, self.unique, self.rescued]
		return "\t".join([str(_) for _ in line])

	def file_header(self, description=""):
		return "\t".join(["#SN", "start", "end", "prob", "unique", "rescued"])

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if parts[0] == "#SN" or len(parts) <= 1:
			return None

		if len(parts) != 6 and len(parts) > 1:
			msg = "Unexpected number of columns in Finesplice junction file.  Expected 6, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[0]
		self.start = int(parts[1])
		self.end = int(parts[2]) - 1

		if fullparse:
			self.score = float(parts[3])
			self.unique = int(parts[4])
			self.rescued = int(parts[5])

		return self


class TruesightJunction(Junction):
	def __init__(self, junc_to_copy=None):
		Junction.__init__(self, use_strand=False, junc_to_copy=junc_to_copy)

		if junc_to_copy and type(junc_to_copy) is TruesightJunction:
			self.mapping = junc_to_copy.mapping
		else:
			self.mapping = 0

	def __str__(self):
		can = 1 if self.canonical == "C" else 2 if self.canonical == "S" else 0
		line = [self.refseq, self.start + 1, self.end + 2, can, self.mapping, self.score]
		return "\t".join([str(_) for _ in line])

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if len(parts) <= 1:
			return None

		if len(parts) != 6 and len(parts) > 1:
			msg = "Unexpected number of columns in Truesight junction file.  Expected 6, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[0]
		self.start = int(parts[1]) - 1
		self.end = int(parts[2]) - 2

		if fullparse:
			self.can = "C" if parts[3] == "1" else "S" if parts[3] == "2" else "N" if parts[3] == "0" else ""
			self.mapping = int(parts[4])
			self.score = float(parts[5])

		return self


class SoapspliceJunction(Junction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		Junction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

	def __str__(self):
		str = "fwd" if self.strand == "+" else "rev" if self.strand == "-" else "fwd"
		line = [self.refseq, self.start, self.end + 2, str, self.score]
		return "\t".join([str(_) for _ in line])

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if len(parts) <= 1:
			return None

		if len(parts) != 5 and len(parts) > 1:
			msg = "Unexpected number of columns in Soapsplice junction file.  Expected 5, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[0]
		self.start = int(parts[1])
		self.end = int(parts[2]) - 2
		self.strand = "+" if parts[3] == "fwd" else "-" if parts[3] == "rev" else "."

		if fullparse:
			self.score = int(parts[4])

		return self


class MapspliceJunction(Junction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		Junction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)

	def __str__(self):
		raise ValueError("Can't represent Mapsplice junctions as string")


	def parse_line(self, line, fullparse=True):
		parts = line.split("\t")

		if len(parts) <= 1:
			return None

		if len(parts) != 29 and len(parts) > 1:
			msg = "Unexpected number of columns in Mapsplice junction file.  Expected 29, found " + str(len(parts))
			raise ValueError(msg)

		self.refseq = parts[0]
		self.start = int(parts[1])
		self.end = int(parts[2]) - 1
		self.strand = parts[5]

		if fullparse:
			self.score = int(parts[4])

		return self


class SpankiJunction(Junction):
	def __init__(self, use_strand=True, junc_to_copy=None):
		Junction.__init__(self, use_strand=use_strand, junc_to_copy=junc_to_copy)
		if junc_to_copy and type(junc_to_copy) is SpankiJunction:
			self.dinucleotide = junc_to_copy.dinucleotide
			self.intron_size = junc_to_copy.intron_size
			self.annotated = junc_to_copy.annotated
			self.rest = junc_to_copy.rest
		else:
			self.dinucleotide = ""
			self.intron_size = 0
			self.annotated = ""
			self.rest = ""

	def __str__(self):
		id_field = self.refseq + ":" + str(self.start) + "_" + str(self.end) + ":" + self.strand
		line = [id_field, self.dinucleotide, self.intron_size, self.annotated, self.rest]
		return "\t".join([str(_) for _ in line])

	def parse_line(self, line, fullparse=True):

		parts = line.split("\t")

		if parts[0] == "juncid" or len(parts) <= 1:
			return None

		if len(parts) != 24 and len(parts) > 1:
			msg = "Unexpected number of columns in Spanki junction file.  Expected 24, found " + str(len(parts))
			raise ValueError(msg)

		parts1 = parts[0].split(":")
		parts2 = parts1[1].split("_")

		self.refseq = parts1[0]
		self.start = int(parts2[0]) - 1
		self.end = int(parts2[1]) - 1
		self.strand = parts1[2]

		if fullparse:
			self.dinucleotide = parts[1]
			self.intron_size = int(parts[2])
			self.annotated = parts[3]
			self.rest = parts[4:]
			self.score = int(parts[9])

		return self


def filterbed(filepath, refset, mode, usestrand, outfile):
	o = open(outfile, 'w')

	o.write("track name=\"junctions\"\n")

	index = 0
	with open(filepath) as i:

		i.readline()  # Skip header
		for line in i:

			key = BedJunction(use_strand=usestrand).parse_line(line, fullparse=False).key

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
