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

#include <intron.hpp>

using portculis::Intron;
using portculis::POSITIVE;
using portculis::NEGATIVE;

BOOST_AUTO_TEST_SUITE(intron)

BOOST_AUTO_TEST_CASE(equality1)
{
    Intron l1(5, 10, 20, POSITIVE);
    Intron l2(5, 10, 20, POSITIVE);
    
    BOOST_CHECK(l1 == l2);
}

BOOST_AUTO_TEST_CASE(equality2)
{
    Intron l1(2, 20, 30, POSITIVE);
    Intron l2(5, 20, 30, POSITIVE);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(equality3)
{
    Intron l1(5, 5, 25, POSITIVE);
    Intron l2(5, 10, 20, POSITIVE);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(equality4)
{
    Intron l1(5, 10, 20, POSITIVE);
    Intron l2(5, 10, 20, NEGATIVE);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(min_anchor) {
    
    Intron intron(5, 10, 20, portculis::POSITIVE);
    
    BOOST_CHECK(intron.minAnchorLength(4, 40) == 6);    
}

BOOST_AUTO_TEST_SUITE_END()
