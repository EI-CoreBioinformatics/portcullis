#!/usr/bin/env python3
__author__ = 'maplesod'


# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

class BedEntry:
	__slots__ = ['__use_strand', 'chrom', 'start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end', 'red',
				 'green', 'blue', 'block_count', 'block_sizes', 'block_starts']

	def __init__(self, use_strand=True):
		self.__use_strand = use_strand
		self.chrom = ""
		self.start = 0
		self.end = 0
		self.name = ""
		self.score = 0
		self.strand = "?"
		self.thick_start = 0
		self.thick_end = 0
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

		line = [self.chrom, self.start, self.end, self.name, self.score,
				self.strand if self.strand else "?",
				self.thick_start, self.thick_end,
				rgb,
				self.block_count,
				bsizes,
				bstarts
				]
		return "\t".join([str(_) for _ in line])

	def __cmp__(self, other):
		if hasattr(other, 'chrom') and hasattr(other, 'start') and hasattr(other, 'end'):
			if self.__lt__(other):
				return 1
			elif self.__gt__(other):
				return -1
			else:
				return 0

	def __key__(self):
		return (self.chrom.encode(), self.thick_start, self.thick_end, self.strand if self.__use_strand else None)

	def __hash__(self):
		return hash(self.__key__())

	@property
	def key(self):
		return self.__key__()

	def __lt__(self, other):
		if self.chrom.__lt__(other.chrom):
			return True
		else:
			if self.chrom.__gt__(other.chrom):
				return False
			else:
				# The same chrom
				if self.thick_start < other.thick_start:
					return True
				else:
					if self.thick_start > other.thick_start:
						return False
					else:
						# Also the same start
						if self.thick_end < other.thick_end:
							return True
						else:
							if self.thick_end > other.thick_end:
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

	@staticmethod
	def create_from_line(key, use_strand=True, tophat=False):

		b = BedEntry(use_strand=use_strand)

		parts = key.split("\t")

		# Handle header or blank lines
		if (len(parts) != 12):
			return None

		b.chrom = parts[0]
		b.start = int(parts[1])
		b.end = int(parts[2])
		b.name = parts[3]
		b.score = int(parts[4])
		b.strand = parts[5]
		b.thick_start = int(parts[6])
		b.thick_end = int(parts[7])

		c_parts = parts[8].split(",")
		b.red = int(c_parts[0])
		b.green = int(c_parts[1])
		b.blue = int(c_parts[2])
		b.block_count = int(parts[9])

		b.block_sizes = [int(_) for _ in parts[10].split(",")]
		# for bp in bsize_parts:
		#     b.block_sizes.append(int(bp))

		b.block_starts = [int(_) for _ in parts[11].rstrip().split(",")]

		if tophat:
			b.thick_start += b.block_sizes[0]
			b.thick_end -= b.block_sizes[1]

		# for bp in bstart_parts:
		#     b.block_starts.append(int(bp))
		assert len(b.block_sizes) == len(b.block_starts) == b.block_count, (key,
																			b.block_count,
																			b.block_sizes,
																			b.block_starts)

		return b

	@staticmethod
	def create_from_tabline(key, use_strand=True, tophat=False):

		b = BedEntry(use_strand=use_strand)

		parts = key.split("\t")

		b.chrom = parts[2]
		b.start = int(parts[6])
		b.end = int(parts[7]) + 1
		b.strand = parts[12]
		b.score = int(parts[14])
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


def makekey(line, usestrand, tophat):
	words = line.split()
	overhang = words[10]
	overhang_parts = overhang.split(",")
	lo = int(overhang_parts[0])
	ro = int(overhang_parts[1])
	chr = words[0]
	start = str(int(words[6]) + lo) if tophat else words[6]
	end = str(int(words[7]) - ro) if tophat else str(int(words[7]))
	strand = words[5]
	if usestrand:
		key = (chr, start, end, strand)
	else:
		key = (chr, start, end, '?')
	return key


def makekeyfromtab(line, usestrand):
	words = line.split("\t")
	chr = words[2]
	start = words[4]
	end = str(int(words[5]) - 1)
	strand = words[11]
	if usestrand:
		key = (chr, start, end, strand)
	else:
		key = (chr, start, end, None)
	return key


def loadbed(filepath, usestrand, tophat):
	index = 0
	items = set()

	with open(filepath) as f:
		# Skip header
		f.readline()
		for line in f:
			key = BedEntry.create_from_line(line,
											use_strand=usestrand,
											tophat=tophat)
			items.add(key)
			index += 1
	if len(items) != index:
		print("duplicated items in bed file " + filepath)
	return items


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
