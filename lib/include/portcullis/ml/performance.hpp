//  ********************************************************************
//  This file is part of Portcullis.
//
//  Portcullis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portcullis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portcullis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::ostream;
using std::stringstream;
using std::shared_ptr;
using std::string;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

namespace portcullis {
namespace ml {

class Performance {
protected:
	uint32_t tp;
	uint32_t tn;
	uint32_t fp;
	uint32_t fn;
public:
	Performance(uint32_t tp, uint32_t tn, uint32_t fp, uint32_t fn) : tp(tp), tn(tn), fp(fp), fn(fn) {}

	inline uint32_t getAllPositive() const {
		return tp + fp;
	}

	inline uint32_t getAllNegative() const {
		return tn + fn;
	}

	inline uint32_t getAllTrue() const {
		return tn + tp;
	}

	inline uint32_t getAllFalse() const {
		return fn + fp;
	}

	inline uint32_t getRealPositive() const {
		return tp + fn;
	}

	inline uint32_t getRealNegative() const {
		return fp + tn;
	}

	inline uint32_t getAll() const {
		return tp + tn + fp + fn;
	}

	inline double getPrecision() const {
		if (getAllPositive() == 0)
			return 0.0;
		else
			return 100.0 * (double)tp / (double)(getAllPositive());
	}

	inline double getPositivePredictiveValue() const {
		return getPrecision();
	}

	inline double getRecall() const {
		if (getRealPositive() == 0)
			return 0.0;
		else
			return 100.0 * (double)tp / (double)(getRealPositive());
	}

	inline double getSensitivity() const {
		return getRecall();
	}

	inline double getTruePositiveRate() const {
		return getRecall();
	}

	inline double getSpecificity() const {
		if (getRealNegative() == 0)
			return 0.0;
		else
			return 100.0 * (double)tn / (double)(getRealNegative());
	}

	inline double getTrueNegativeRate() const {
		return getSpecificity();
	}

	inline double getNPV() const {
		if (getAllNegative() == 0)
			return 0.0;
		else
			return 100.0 * (double)tn / (double)(getAllNegative());
	}

	inline double getFallOut() const {
		return 100.0 - getSpecificity();
	}

	inline double getFalsePositiveRate() const {
		return getFallOut();
	}

	inline double getFDR() const {
		return 100.0 - getPrecision();
	}

	inline double getFNR() const {
		return 100.0 * getRecall();
	}

	/**
	 * The proportion of the population that really have the positive condition
	 * (TP+FN).  100% means the whole population has the positive condition.  50%
	 * means half do and 0% means none do.  This value is not effected by the predictor.
	 */
	inline double getPrevalence() const {
		if (getAll() == 0)
			return 0.0;
		else
			return 100.0 * (double)getRealPositive() / (double)getAll();
	}

	/**
	 * The bias between positive predictions (TP + FP) and the population.  This
	 * represents the bias of the system towards positive predictions and can
	 * therefore be modified by adjusting the model.
	 */
	inline double getBias() const {
		if (getAll() == 0)
			return 0.0;
		else
			return 100.0 * (double)getAllPositive() / (double)getAll();
	}

	/**
	 * The overall accuracy of the system, representing the closeness to the true
	 * sample.  Note however that this is a biased metric in that it does not
	 * properly cater for prevalence of positives conditions and bias of the model.
	 */
	inline double getAccuracy() const {
		if (getAll() == 0)
			return 0.0;
		else
			return 100.0 * (double)getAllTrue() / (double)getAll();
	}

	/**
	 * Returns the F-score with provided beta parameter.  F-Scores are biased in
	 * favour of positive values: true negatives can vary without impacting the
	 * F-Score.
	 * @param beta A positive real number.  1.0 sets a balance between precision
	 * and recall.  < 1.0 biases in favour or precison.  > 1.0 biases in favour
	 * of recall.
	 */
	inline double getFBScore(const double beta) const {
		if (beta <= 0) return 0.0;
		const double recall = getRecall();
		const double precision = getPrecision();
		const double beta2 = beta * beta;
		return (double)(1.0 + beta2) * (precision * recall) / ((beta2 * precision) + recall);
	}

	/**
	 * Returns the balanced F-score, which is the harmonic mean of precision and
	 * recall.
	 */
	inline double getF1Score() const {
		return getFBScore(1.0);
	}

	/**
	 * Matthew's Correlation Coefficient.  An unbiased measure generally considered
	 * one of the best measures of overall system performance.
	 */
	inline double getMCC() const {
		return std::sqrt(getInformedness() * getMarkedness());
	}

	/**
	 * Informedness is an unbiased measure that gives the probability that
	 * you have made an informed decision (versus chance).
	 */
	inline double getInformedness() const {
		return getSensitivity() + getSpecificity() - 100.0;
	}

	/**
	 * Markedness is an unbiased measure that gives the probability that a
	 * condition is marked by the predictor (versus chance).
	 */
	inline double getMarkedness() const {
		return getPrecision() + getNPV() - 100.0;
	}

	string toShortString() const;
	string toLongString() const;

	static string shortHeader() {
		return "TP\tTN\tFP\tFN\tREC\tPRC\tF1";
	}

	static string longHeader() {
		return "TP\tTN\tFP\tFN\tPREV\tBIAS\tSENS\tSPEC\tPPV\tNPV\tF1\tACC\tINFO\tMARK\tMCC";
	}

	static string to_2dp_string(const double v) {
		std::ostringstream out;
		out << std::fixed << std::setprecision(2) << v;
		return out.str();
	}

	static void loadGenuine(path& genuineFile, vector<bool>& results);

};

class PerformanceList {
public:

	void clear() {
		this->scores.clear();
	}

	shared_ptr<Performance> operator [](int i) const {
		return this->scores[i];
	}

	void add(shared_ptr<Performance> p) {
		this->scores.push_back(p);
	}

	void outputMeanPerformance(std::ostream& resout);

protected:

	vector<shared_ptr<Performance>> scores;

	void outputMeanScore(const vector<double>& scores, const string& score_type, std::ostream& resout);

};

}
}