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

#include <iostream>
using std::cout;
using std::endl;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>

#include <portcullis/seq_utils.hpp>
using portcullis::SeqUtils;

#include <portcullis/ml/markov_model.hpp>
using portcullis::ml::KMMU;

void portcullis::ml::KmerMarkovModel::train(const vector<string>& input, const uint32_t _order) {
	order = _order;
	KMMU temp;
	for (auto & seq : input) {
		string s = SeqUtils::makeClean(seq);
		if (s.size() > order + 1) {
			for (size_t i = order; i < s.size(); i++) {
				temp[s.substr(i - order, order)][s.substr(i, 1)]++;
			}
		}
	}
	model.clear();
	for (auto & i : temp) {
		double sum = 0;
		for (auto & j : i.second) {
			sum += j.second;
		}
		for (auto & j : i.second) {
			//cout << i.first << " " << j.first << " " << j.second << endl;
			model[i.first][j.first] = j.second / sum;
		}
	}
}


double portcullis::ml::KmerMarkovModel::getScore(const string& seq) {
	string s = SeqUtils::makeClean(seq);
	double score = 1.0;
	uint32_t no_count = 0;
	for (size_t i = order; i < s.size(); i++) {
		double m = model[s.substr(i - order, order)][s.substr(i, 1)];
		if (m != 0.0) {
			score *= m;
		}
		else {
			no_count++;
		}
	}
	if (score == 0.0) {
		return -100.0;
	}
	else if (no_count > 2) {
		// Add a penalty for situations where we repeatedly don't find a kmer in the tranining set
		score /= ((double)no_count * 0.5);
	}
	return log(score);
}

void portcullis::ml::PosMarkovModel::train(const vector<string>& input, const uint32_t _order) {
	order = _order;
	PMMU temp;
	for (auto & seq : input) {
		string s = SeqUtils::makeClean(seq);
		for (size_t i = order; i < s.size(); i++) {
			temp[i][s.substr(i, 1)]++;
		}
	}
	model.clear();
	for (auto & i : temp) {
		double sum = 0;
		for (auto & j : i.second) {
			sum += j.second;
		}
		for (auto & j : i.second) {
			//cout << i.first << " " << j.first << " " << j.second << endl;
			model[i.first][j.first] = j.second / sum;
		}
	}
}


double portcullis::ml::PosMarkovModel::getScore(const string& seq) {
	string s = SeqUtils::makeClean(seq);
	double score = 1.0;
	for (size_t i = order; i < s.size(); i++) {
		score *= model[i][s.substr(i, 1)];
	}
	if (score == 0.0) {
		return -300.0;
	}
	return log(score);
}

