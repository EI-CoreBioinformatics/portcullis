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

#include <location.hpp>

using std::cout;
using std::endl;

using portculis::Location;

BOOST_AUTO_TEST_SUITE(location)

BOOST_AUTO_TEST_CASE(equality1)
{
    Location l1(5, 10, 20);
    Location l2(5, 10, 20);
    
    BOOST_CHECK(l1 == l2);
}

BOOST_AUTO_TEST_CASE(equality2)
{
    Location l1(2, 20, 30);
    Location l2(5, 20, 30);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_CASE(equality3)
{
    Location l1(5, 5, 25);
    Location l2(5, 10, 20);
    
    BOOST_CHECK(l1 != l2);
}

BOOST_AUTO_TEST_SUITE_END()