#!/usr/bin/env python3
__author__ = 'maplesod'


class PEntry:
	def __init__(self):
		self.tp = 0
		self.fp = 0
		self.fn = 0

	def __str__(self):
		return str(self.tp) + "\t" + str(self.fp) + "\t" + str(self.fn) + "\t" + format(self.calcRecall(), '.2f') \
			   + "\t" + format(self.calcPrecision(), '.2f') + "\t" + format(self.calcF1(), '.2f')

	def calcPrecision(self):
		p = float(self.tp + self.fp)
		return 100.0 * float(self.tp) / p if p > 0.0 else 0.0

	def calcRecall(self):
		tpfn = float(self.tp + self.fn)
		return 100.0 * float(self.tp) / tpfn if tpfn > 0.0 else 0.0

	def calcF1(self):
		prc = self.calcPrecision()
		rec = self.calcRecall()
		return 2.0 * (prc * rec) / (prc + rec) if (prc + rec) > 0.0 else 0.0

	@staticmethod
	def header():
		return "TP\tFP\tFN\tRecall\tPrecision\tF1"
