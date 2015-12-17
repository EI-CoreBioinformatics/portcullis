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

#include "../src/intron.hpp"
#include "../src/junction.hpp"
using portcullis::CanonicalSS;
using portcullis::Intron;
using portcullis::Junction;
using portcullis::JunctionException;

bool is_critical( JunctionException const& ex ) { return true; }

const RefSeq rd2(2, "seq_2", 100);
const RefSeq rd5(5, "seq_5", 100);


TEST(junction, intron) {
    
    shared_ptr<Intron> l1(new Intron(rd5, 20, 30));
    Junction j1(l1, 10, 40);
    
    int32_t intronSz = j1.getIntronSize();
    EXPECT_EQ(intronSz, 11);
}

TEST(junction, donor_acceptor) {
    
    shared_ptr<Intron> l1(new Intron(rd5, 20, 30));
    Junction j1(l1, 10, 40);
    
    shared_ptr<Intron> l2(new Intron(rd5, 20, 30));
    Junction j2(l2, 10, 40);
    
    CanonicalSS res1 = j1.setDonorAndAcceptorMotif("GT", "AG");
    EXPECT_EQ(res1, portcullis::CANONICAL);
    
    CanonicalSS res2 = j2.setDonorAndAcceptorMotif("CT", "AC");
    EXPECT_EQ(res2, portcullis::CANONICAL);
    
    try {
        j1.setDonorAndAcceptorMotif("GTA", "AG");
    }
    catch(JunctionException& err) {
        EXPECT_EQ(true, true);
    }
    catch(...) {
        FAIL() << "Expected Junction Exception";
    }

    
    CanonicalSS res4 = j1.setDonorAndAcceptorMotif("CT", "AG");
    EXPECT_NE(res4, portcullis::CANONICAL);
    
    CanonicalSS res5 = j1.setDonorAndAcceptorMotif("GT", "AC");
    EXPECT_NE(res5, portcullis::CANONICAL);
    
    try {
        j1.setDonorAndAcceptorMotif("", "");
    }
    catch(JunctionException& err) {
        EXPECT_EQ(true, true);
    }
    catch(...) {
        FAIL() << "Expected Junction Exception";
    }
}

TEST(junction, entropy) {
    
    shared_ptr<Intron> l(new Intron(rd5, 20, 30));
    Junction j(l, 10, 40);
    
    int32_t ints1[] = {13, 15, 17, 19};
    vector<int32_t> juncPos1(ints1, ints1 + sizeof(ints1) / sizeof(int32_t)); 
    
    int32_t ints2[] = {16, 16, 16, 16};
    vector<int32_t> juncPos2(ints2, ints2 + sizeof(ints2) / sizeof(int32_t)); 
    
    double e1 = j.calcEntropy(juncPos1);
    double e2 = j.calcEntropy(juncPos2);
    
    // Actually I don't know what the entropy scores should be exactly... but e1 
    // should definitely have a higher entropy than e2
    EXPECT_GT(e1, e2);
}

/**
 * This IS what you'd expect to see in a real junction
 */
TEST(junction, coverage1) {
    
    shared_ptr<Intron> l(new Intron(rd5, 20, 30));
    Junction j1(l, 10, 40);
    
    vector<uint32_t> coverage1{ 10,10,10,10,10,10,10,10,10,10,
                                10,10,10,10,10,8,6,4,3,2,
                                0,0,0,0,0,0,0,0,0,0,
                                2,3,4,7,8,10,10,10,10,10,
                                10,10,10,10,10,10,10,10,10,10}; 
    
    double cvg1 = j1.calcCoverage(coverage1);
                               
    //cout << "Coverage: " << cvg1 << endl;
    EXPECT_GT(cvg1, 0);
}

/**
 * This IS NOT what you'd expect to see in a real junction
 */
TEST(junction, coverage2) {
    
    shared_ptr<Intron> l(new Intron(rd5, 20, 30));
    Junction j2(l, 10, 40);    
    
    vector<uint32_t> coverage2{ 0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,2,3,5,7,8,
                                10,10,10,10,10,10,10,10,10,10,
                                8,6,4,3,2,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0}; 
    
    double cvg2 = j2.calcCoverage(coverage2);
                               
    //cout << "Coverage: " << cvg2 << endl;
    
    EXPECT_LT(cvg2, 0);
}
