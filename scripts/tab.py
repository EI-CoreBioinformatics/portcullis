#!/usr/bin/env python3

import bed12

class TabEntry:
	chrom = ""
	start = 0
	end = 0
	left = 0
	right = 0
	strand = "?"
	#M1 = "N"
	M2 = 0
	M3 = 0
	M4 = 0
	#M5 = 0
	#M6 = 0
	#M7 = 0
	M8 = 0
	M9 = 0
	M10 = 0
	M11 = 0.0
	M12 = 0
	M13 = 0
	M14 = 0
	#M15 = 0.0

	def __init__(self):
		self.data = []

	def __str__(self):

		line = [ self.chrom, self.start, self.end, self.left, self.right, self.strand,
			 self.M2, self.M3, self.M4, self.M8, self.M9, self.M10, self.M11, self.M12, self.M13, self.M14
			 ]
		return "\t".join([str(_) for _ in line])

	def makeMatrixRow(self):
		return [ self.M2, self.M3, self.M4, self.M8, self.M9, self.M10, self.M11, self.M12, self.M13, self.M14 ]

	@staticmethod
	def features():
		return [ "M2-nb-reads", "M3-nb_dist_aln", "M4-nb_rel_aln", "M8-max_min_anc", "M9-dif_anc", "M10-dist_anc",
				 "M11-entropy", "M12-maxmmes", "M13-hammping5p", "M14-hamming3p"]
	@staticmethod
	def featureAt(index):
		return TabEntry.features()[index]

	@staticmethod
	def sortedFeatures(indicies):
		f = TabEntry.features()
		s=[]
		for i in indicies:
			s.append(f[i])
		return s

	@staticmethod
	def nbMetrics():
		return 10

	@staticmethod
	def create_from_tabline(line):

		b = TabEntry()

		parts = line.split("\t")

		b.chrom = parts[2]
		b.left = int(parts[6])
		b.right = int(parts[7])
		b.strand = parts[12]
		b.start = int(parts[4])
		b.end = int(parts[5])

		b.M2 = int(parts[14])
		b.M3 = int(parts[15])
		b.M4 = int(parts[16])
		b.M8 = int(parts[20])
		b.M9 = int(parts[21])
		b.M10 = int(parts[22])
		b.M11 = float(parts[23])
		b.M12 = int(parts[24])
		b.M13 = int(parts[25])
		b.M14 = int(parts[26])

		return b



def loadtab(tabfile):

	bed=list()
	tab=list()
	with open(tabfile) as f:
		# Skip header
		f.readline()

		for line in f:
			line.strip()
			if len(line) > 1:
				bed.append(bed12.BedEntry.create_from_tabline(line, False, False))
				tab.append(TabEntry.create_from_tabline(line))


	return bed, tab

