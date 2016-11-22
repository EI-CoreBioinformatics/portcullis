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

#include <algorithm>
#include <iterator>
#include <vector>
using std::vector;
using std::random_shuffle;

#include <boost/exception/all.hpp>

namespace portcullis {
namespace ml {

typedef boost::error_info<struct KFoldError, string> KFoldErrorInfo;
struct KFoldException: virtual boost::exception, virtual std::exception { };

// Derived from https://sureshamrita.wordpress.com/2011/08/24/c-implementation-of-k-fold-cross-validation/
template<class In>
class KFold {
public:
	KFold(int k, In _beg, In _end) :
		beg(_beg), end(_end), K(k) {
		if (K <= 0)
			BOOST_THROW_EXCEPTION(KFoldException() << KFoldErrorInfo(string(
									  "The supplied value of K is =") + lexical_cast<string>(K) +
								  ". One cannot create " + lexical_cast<string>(K) + "no of folds"));
		//create the vector of integers
		int foldNo = 0;
		for (In i = beg; i != end; i++) {
			whichFoldToGo.push_back(++foldNo);
			if (foldNo == K)
				foldNo = 0;
		}
		if (!K)
			BOOST_THROW_EXCEPTION(KFoldException() << KFoldErrorInfo(string(
									  "With this value of k (=") + lexical_cast<string>(K) +
								  ")Equal division of the data is not possible"));
		random_shuffle(whichFoldToGo.begin(), whichFoldToGo.end());
	}

	template<class Out>
	void getFold(int foldNo, Out training, Out testing) {
		int k = 0;
		In i = beg;
		while (i != end) {
			if (whichFoldToGo[k++] == foldNo) {
				*testing++ = *i++;
			}
			else
				*training++ = *i++;
		}
	}

private:
	In beg;
	In end;
	int K; //how many folds in this
	vector<int> whichFoldToGo;
};

}
}