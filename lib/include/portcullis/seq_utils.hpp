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

#include <string>
using std::string;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;


namespace portcullis {
    
const char REVCOMP_LOOKUP[] = {'T',  0,  'G', 'H',
                                0,   0,  'C', 'D',
                                0,   0,   0,   0,
                               'K', 'N',  0,   0,
                                0,  'Y', 'W', 'A',
                               'A', 'B', 'S', 'X',
                               'R',  0 };

typedef boost::error_info<struct SeqUtilsError,string> SeqUtilsErrorInfo;
struct SeqUtilsException: virtual boost::exception, virtual std::exception { };




class SeqUtils {
    
public:
    
    static int16_t hammingDistance(const string& s1, const string& s2) {

        if (s1.size() != s2.size())
            BOOST_THROW_EXCEPTION(SeqUtilsException() << SeqUtilsErrorInfo(string(
                    "Can't find hamming distance of strings that are not the same length.  ") +
                    "s1: " + lexical_cast<string>(s1.size()) + "\"" + s1 + "\"; " +
                    "s2: " + lexical_cast<string>(s2.size()) + "\"" + s2 + "\""));

        string s1u = boost::to_upper_copy(s1);
        string s2u = boost::to_upper_copy(s2);

        int16_t sum = 0;
        for(size_t i = 0; i < s1u.size(); i++) {
            if (s1u[i] != s2u[i]) {
                sum++;
            }
        }

        return sum;
    }

    /**
     * Reverses a string, and returns a copy of that string
     * @param sequence
     * @return 
     */
    static string reverseSeq(string& sequence) {
        return string(sequence.rbegin(), sequence.rend());
    }

    /**
     * Returns a reverse complement of the provided sequence
     * @param sequence
     * @return 
     */
    static string reverseComplement(string sequence) {

        // do complement, in-place
        size_t seqLength = sequence.length();
        for ( size_t i = 0; i < seqLength; ++i )
            sequence.replace(i, 1, 1, REVCOMP_LOOKUP[(int)sequence.at(i) - 65]);

        // reverse it
        return reverseSeq(sequence);
    }
};
    
    
}

