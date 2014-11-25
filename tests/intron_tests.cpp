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
#include <boost/test/unit_test.hpp>

#include <intron.hpp>
using portcullis::RefSeq;
using portcullis::Intron;
using portcullis::POSITIVE;
using portcullis::NEGATIVE;

const RefSeq rd2(2, "seq_2", 100);
const RefSeq rd5(5, "seq_5", 100);

BOOST_AUTO_TEST_SUITE(intron)

BOOST_AUTO_TEST_CASE(equality1)
{
    Intron l1(rd5, 10, 20, POSITIVE);
    Intron l2(rd5, 10, 20, POSITIVE);
    
    BOOST_CHECK(l1 == l2);
}

BOOST_AUTO_TEST_CASE(equality2)
{
    Intron l1(rd2, 20, 30, POSITIVE);
    Intron l2(rd5, 20, 30, POSITIVE);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(equality3)
{
    Intron l1(rd5, 5, 25, POSITIVE);
    Intron l2(rd5, 10, 20, POSITIVE);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(equality4)
{
    Intron l1(rd5, 10, 20, POSITIVE);
    Intron l2(rd5, 10, 20, NEGATIVE);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(min_anchor) {
    
    Intron intron(rd5, 10, 20, portcullis::POSITIVE);
    
    BOOST_CHECK(intron.minAnchorLength(4, 40) == 6);    
}

BOOST_AUTO_TEST_CASE(size) {
    
    Intron intron(rd5, 10, 20, portcullis::POSITIVE);
    
    BOOST_CHECK(intron.size() == 11);    
}

BOOST_AUTO_TEST_SUITE_END()
