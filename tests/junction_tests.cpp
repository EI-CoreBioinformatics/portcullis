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

using portculis::Location;
using portculis::Junction;

BOOST_AUTO_TEST_SUITE(junction)

BOOST_AUTO_TEST_CASE(intron)
{
    shared_ptr<Location> l1(new Location(5, 20, 30));
    Junction j1(l1, 10, 40);
    
    int32_t intronSz = j1.getIntronSize();
    BOOST_CHECK(intronSz == 10);
}

BOOST_AUTO_TEST_CASE(donor_acceptor)
{
    shared_ptr<Location> l1(new Location(5, 20, 30));
    Junction j1(l1, 10, 40);
    
    bool res1 = j1.setDonorAndAcceptorMotif("GT", "AG");
    BOOST_CHECK(res1);
    
    bool res2 = j1.setDonorAndAcceptorMotif("CT", "AC");
    BOOST_CHECK(res2);
    
    bool res3 = j1.setDonorAndAcceptorMotif("GTA", "AG");
    BOOST_CHECK(!res3);
    
    bool res4 = j1.setDonorAndAcceptorMotif("CT", "AG");
    BOOST_CHECK(!res4);
    
    bool res5 = j1.setDonorAndAcceptorMotif("GT", "AC");
    BOOST_CHECK(!res5);
    
    bool res6 = j1.setDonorAndAcceptorMotif("", "");
    BOOST_CHECK(!res6);        
}

BOOST_AUTO_TEST_SUITE_END()