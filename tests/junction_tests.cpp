//  ********************************************************************
//  This file is part of Portculis.
//
//  Portculis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portculis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE PORTCULIS
#endif
#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>

#include <junction.hpp>

using std::cout;
using std::endl;

using portculis::Intron;
using portculis::Junction;
using portculis::JunctionException;

bool is_critical( JunctionException const& ex ) { return true; }

BOOST_AUTO_TEST_SUITE(junction)

BOOST_AUTO_TEST_CASE(intron) {
    
    shared_ptr<Intron> l1(new Intron(5, 20, 30, portculis::POSITIVE));
    Junction j1(l1, 10, 40);
    
    int32_t intronSz = j1.getIntronSize();
    BOOST_CHECK(intronSz == 11);
}

BOOST_AUTO_TEST_CASE(donor_acceptor) {
    
    shared_ptr<Intron> l1(new Intron(5, 20, 30, portculis::POSITIVE));
    Junction j1(l1, 10, 40);
    
    shared_ptr<Intron> l2(new Intron(5, 20, 30, portculis::NEGATIVE));
    Junction j2(l2, 10, 40);
    
    bool res1 = j1.setDonorAndAcceptorMotif("GT", "AG");
    BOOST_CHECK(res1);
    
    bool res2 = j2.setDonorAndAcceptorMotif("CT", "AC");
    BOOST_CHECK(res2);
    
    BOOST_CHECK_EXCEPTION(j1.setDonorAndAcceptorMotif("GTA", "AG"), JunctionException, is_critical);
    
    bool res4 = j1.setDonorAndAcceptorMotif("CT", "AG");
    BOOST_CHECK(!res4);
    
    bool res5 = j1.setDonorAndAcceptorMotif("GT", "AC");
    BOOST_CHECK(!res5);
    
    BOOST_CHECK_EXCEPTION(j1.setDonorAndAcceptorMotif("", ""), JunctionException, is_critical);
}

BOOST_AUTO_TEST_CASE(entropy) {
    
    shared_ptr<Intron> l(new Intron(5, 20, 30, portculis::POSITIVE));
    Junction j(l, 10, 40);
    
    int32_t ints1[] = {13, 15, 17, 19};
    vector<int32_t> juncPos1(ints1, ints1 + sizeof(ints1) / sizeof(int32_t)); 
    
    int32_t ints2[] = {16, 16, 16, 16};
    vector<int32_t> juncPos2(ints2, ints2 + sizeof(ints2) / sizeof(int32_t)); 
    
    double e1 = j.calcEntropy(juncPos1);
    double e2 = j.calcEntropy(juncPos2);
    
    // Actually I don't know what the entropy scores should be exactly... but e1 
    // should definitely have a higher entropy than e2
    BOOST_CHECK(e1 > e2);
}

BOOST_AUTO_TEST_SUITE_END()