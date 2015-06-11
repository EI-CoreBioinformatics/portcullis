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

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE PORTCULLIS
#endif

#include <iostream>
using std::cout;
using std::endl;

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "seq_utils.hpp"
using portcullis::SeqUtils;

BOOST_AUTO_TEST_SUITE(seq_utils)

BOOST_AUTO_TEST_CASE(hamming) {
    
    BOOST_CHECK(SeqUtils::hammingDistance("ATGC", "ATGC") == 0);
    BOOST_CHECK(SeqUtils::hammingDistance("ATGC", "ATGG") == 1);
    BOOST_CHECK(SeqUtils::hammingDistance("ATGC", "CGTA") == 4);
}


BOOST_AUTO_TEST_CASE(rev) {
    
    string seq("ATGC");
    BOOST_CHECK(SeqUtils::reverseSeq(seq) == "CGTA");
}

BOOST_AUTO_TEST_CASE(rev_comp) {
    
    BOOST_CHECK(SeqUtils::reverseComplement("ATGC") == "GCAT");    
}


BOOST_AUTO_TEST_SUITE_END()
