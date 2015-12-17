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

#include <gtest/gtest.h>

#include <iostream>
using std::cout;
using std::endl;

#include <boost/filesystem.hpp>

#include "../src/seq_utils.hpp"
using portcullis::SeqUtils;


TEST(seq_utils, hamming) {
    
    EXPECT_EQ(SeqUtils::hammingDistance("ATGC", "ATGC"), 0);
    EXPECT_EQ(SeqUtils::hammingDistance("ATGC", "ATGG"), 1);
    EXPECT_EQ(SeqUtils::hammingDistance("ATGC", "CGTA"), 4);
}


TEST(seq_utils, rev) {
    
    string seq("ATGC");
    EXPECT_EQ(SeqUtils::reverseSeq(seq), "CGTA");
}

TEST(seq_utils, rev_comp) {
    
    EXPECT_EQ(SeqUtils::reverseComplement("ATGC"), "GCAT");    
}

