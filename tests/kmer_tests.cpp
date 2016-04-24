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
#include <fstream>
#include <string>
using std::cout;
using std::endl;
using std::stringstream;


#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <portcullis/kmer.hpp>
using portcullis::KmerHash;

        
TEST(kmer, make) {
    
    string seq = "ATGCATGCATCGNATATATATTGAC";
    
    KmerHash kh(4, seq);
    
    //kh.print(std::cout);
    
    EXPECT_EQ(1, kh.getCount("tgac"));
    
    EXPECT_EQ(0, kh.getCount("kjls"));
    
    EXPECT_EQ(16, kh.nbDistinctKmers());
    
}

