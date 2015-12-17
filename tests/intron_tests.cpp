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

#include "../src/bam/bam_master.hpp"
using portcullis::bam::RefSeq;

#include "../src/intron.hpp"
using portcullis::Intron;

const RefSeq rd2(2, "seq_2", 100);
const RefSeq rd5(5, "seq_5", 100);

TEST(intron, equality1) {
    
    Intron l1(rd5, 10, 20);
    Intron l2(rd5, 10, 20);
    
    EXPECT_EQ(l1, l2);
}

TEST(intron, equality2) {
    
    Intron l1(rd2, 20, 30);
    Intron l2(rd5, 20, 30);
    
    EXPECT_NE(l1, l2);
}

TEST(intron, equality3) {
    
    Intron l1(rd5, 5, 25);
    Intron l2(rd5, 10, 20);
    
    EXPECT_NE(l1, l2);
}

TEST(intron, min_anchor) {
    
    Intron intron(rd5, 10, 20);
    
    EXPECT_EQ(intron.minAnchorLength(4, 40), 6);    
}

TEST(intron, size) {
    
    Intron intron(rd5, 10, 20);
    
    EXPECT_EQ(intron.size(), 11);    
}
