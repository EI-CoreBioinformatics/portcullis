#!/usr/bin/env python3

import abc

__author__ = 'maplesod'


class Junction( object ):
	__metaclass__  = abc.ABCMeta

	def __init__(self, use_strand=True):
		self.__use_strand = use_strand
		self.refseq = ""
		self.start = 0
		self.end = 0
		self.strand = "."

	def __str__(self):
		line = [self.refseq, self.start, self.end, self.strand]
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
		return (self.refseq.encode(), self.start, self.end, self.strand if self.__use_strand else None)

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
	def parse_line(line): pass

	def createDict(self, filepath):

		index = 0
		items = set()

		with open(filepath) as f:
			# Skip header
			f.readline()
			for line in f:
				key = self.parse_line(line).key
				items[key] = line
				index += 1
		if len(items) != index:
			print("duplicated items in bed file " + filepath)
		return items


class BedEntry(Junction):

	def __init__(self, use_strand=True, tophat=False):
		Junction.__init__(use_strand)
		self.tophat = tophat
		self.name = ""
		self.score = 0
		self.left = 0
		self.right = 0
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


	def parse_line(self, line):

		parts = line.split("\t")

		# Handle header or blank lines
		if (len(parts) != 12):
			return None

		self.refseq = parts[0]
		self.left = int(parts[1])
		self.right = int(parts[2])
		self.name = parts[3]
		self.score = float(parts[4])
		self.strand = parts[5]
		self.start = int(parts[6])
		self.end = int(parts[7])

		c_parts = parts[8].split(",")
		self.red = int(c_parts[0])
		self.green = int(c_parts[1])
		self.blue = int(c_parts[2])
		self.block_count = int(parts[9])

		self.block_sizes = [int(_) for _ in parts[10].split(",")]
		# for bp in bsize_parts:
		#     b.block_sizes.append(int(bp))

		self.block_starts = [int(_) for _ in parts[11].rstrip().split(",")]

		if self.tophat:
			self.start += self.block_sizes[0]
			self.end -= self.block_sizes[1]

		# for bp in bstart_parts:
		#     b.block_starts.append(int(bp))
		assert len(self.block_sizes) == len(self.block_starts) == self.block_count, (line,
																			self.block_count,
																			self.block_sizes,
																			self.block_starts)

		return self

	@staticmethod
	def create_from_tabline(key, use_strand=True):

		b = BedEntry(use_strand=use_strand)

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

			key = makekey(line, usestrand, tophat)

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
				key = makekeyfromtab(line, usestrand)
				if mode == 0:
					if key in bed:
						o.write(line + "\n")
				elif mode == 1:
					if key not in bed:
						o.write(line + "\n")

	o.close()
