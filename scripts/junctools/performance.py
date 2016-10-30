#!/usr/bin/env python3
__author__ = 'maplesod'

import math

class Performance:
	def __init__(self, tp=0, fp=0, fn=0, tn=0):
		self.tp = tp
		self.fp = fp
		self.fn = fn
		self.tn = tn

	def positives(self):
		return self.tp + self.fp

	def negatives(self):
		return self.tn + self.fn

	def trues(self):
		return self.tn + self.tp

	def falses(self):
		return self.fn + self.fp

	def RP(self):
		return self.tp + self.fn

	def RN(self):
		return self.fp + self.tn

	def N(self):
		return self.tp + self.tn + self.fp + self.fn

	def prevalence(self):
		"""
		The proportion of the population that really have the positive condition
		(TP+FN).  100% means the whole population has the positive condition.  50%
		means half do and 0% means none do.  This value is not effected by the predictor.
		"""
		return 100.0 * float(self.RP()) / float(self.N())

	def bias(self):
		"""
		The bias between positive predictions (TP + FP) and the population.  This
		represents the bias of the system towards positive predictions and can
		therefore be modified by adjusting the model.
		"""
		return 100.0 * float(self.positives()) / float(self.N())

	def precision(self):
		p = float(self.positives())
		return 100.0 * float(self.tp) / p if p > 0.0 else 0.0

	def npv(self):
		n = float(self.negatives())
		return 100.0 * float(self.tn) / n if n > 0.0 else 0.0

	def recall(self):
		rp = float(self.RP())
		return 100.0 * float(self.tp) / rp if rp > 0.0 else 0.0

	def specificity(self):
		rn = float(self.RN())
		return 100.0 * float(self.tn) / rn if rn > 0.0 else 0.0

	def F1(self):
		prc = self.precision()
		rec = self.recall()
		return 2.0 * (prc * rec) / (prc + rec) if (prc + rec) > 0.0 else 0.0

	def accuracy(self):
		"""
		The overall accuracy of the system, representing the closeness to the true sample.  Note however that this is a
		biased metric in that it does not properly cater for prevalence of positives conditions and bias of the model.
		"""
		n = float(self.N())
		return 100.0 * float(self.trues()) / n if n > 0.0 else 0.0

	def informedness(self):
		"""
		Informedness is an unbiased measure that gives the probability that you have made an informed decision (versus chance).
		:param self:
		:return:
		"""
		return self.recall() + self.specificity() - 100.0

	def markedness(self):
		"""
		Markedness is an unbiased measure that gives the probability that a condition is marked by the predictor (versus
		chance).
		:param self:
		:return:
		"""
		return self.precision() + self.npv() - 100.0;

	def MCC(self):
		"""
		Matthew's Correlation Coefficient.  An unbiased measure generally considered one of the best measures of overall
		system performance.
		:param self:
		:return:
		"""
		return math.sqrt(self.informedness() * self.markedness()) if not self.informedness() == 0.0 and not self.markedness() == 0.0 else 0.0


	def __str__(self):
		return self.shortStr()

	def shortStr(self):
		parts = []
		parts.append(str(self.tp))
		parts.append(str(self.fp))
		parts.append(str(self.fn))
		parts.append(format(self.recall(), '.2f'))
		parts.append(format(self.precision(), '.2f'))
		parts.append(format(self.F1(), '.2f'))
		return "\t".join(parts)


	def longStr(self):
		parts = []
		parts.append(str(self.tp))
		parts.append(str(self.tn))
		parts.append(str(self.fp))
		parts.append(str(self.fn))
		parts.append(format(self.prevalence(), '.2f'))
		parts.append(format(self.bias(), '.2f'))
		parts.append(format(self.recall(), '.2f'))
		parts.append(format(self.specificity(), '.2f'))
		parts.append(format(self.precision(), '.2f'))
		parts.append(format(self.npv(), '.2f'))
		parts.append(format(self.F1(), '.2f'))
		parts.append(format(self.accuracy(), '.2f'))
		parts.append(format(self.informedness(), '.2f'))
		parts.append(format(self.markedness(), '.2f'))
		parts.append(format(self.MCC(), '.2f'))

		return "\t".join(parts)


	@staticmethod
	def shortHeader():
		return "TP\tFP\tFN\tREC\tPRC\tF1"


	@staticmethod
	def longHeader():
		return "TP\tTN\tFP\tFN\tPREV\tBIAS\tSENS\tSPEC\tPPV\tNPV\tF1\tACC\tINFO\tMARK\tMCC"
