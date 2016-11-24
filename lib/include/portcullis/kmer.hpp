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

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using std::ostream;
using std::string;
using std::unordered_map;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>

namespace portcullis {

typedef boost::error_info<struct KmerError, string> KmerErrorInfo;
struct KmerException: virtual boost::exception, virtual std::exception { };

/**
 * An extremely basic kmer hash class that doesn't worry about canonical kmers
 * or reverse complementing the data.
 */
class KmerHash {
private:
	unordered_map<string, uint32_t> khash;
	uint16_t k;

public:
	KmerHash(const uint16_t _k, const string& seq) {
		k = _k;
		for (size_t i = k; i <= seq.size(); i++) {
			khash[boost::to_upper_copy(seq.substr(i - k, k))]++;
		}
	}

	uint32_t getCount(const string& kmer) {
		if (kmer.size() != k) {
			BOOST_THROW_EXCEPTION(KmerException() << KmerErrorInfo(string(
									  "Given kmer size is ") + std::to_string(kmer.size()) + ".  Expected kmer of size " + std::to_string(k)));
		}
		string kup = boost::to_upper_copy(kmer);
		return (khash.find(kup) != khash.end()) ? khash[kup] : 0;
	}

	size_t nbDistinctKmers() const {
		return khash.size();
	}

	void print(ostream& out) const {
		for (auto & kmer : khash) {
			out << kmer.first << "\t" << kmer.second << endl;
		}
	}

	void printAbundanceHistogram(ostream& out, const uint32_t hist_size) {
		vector<uint32_t> hist(hist_size, 0);
		for (auto & kmer : khash) {
			if (hist.size() > kmer.second) {
				hist[kmer.second]++;
			}
			else {
				hist[hist.size() - 1]++;
			}
		}
		for (auto & e : hist) {
			out << e << endl;
		}
	}
};

}