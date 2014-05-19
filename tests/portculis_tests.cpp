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

#include <portculis.hpp>

BOOST_AUTO_TEST_SUITE(portculis)

BOOST_AUTO_TEST_CASE(constructor)
{
    portculis::Portculis portculis(
                        "tests/resources/ecoli.bam",
                        "tests/resources/ecoli.fa",
                        "tests/tmp/portculis/constructor1",
                        1,
                        false,
                        false);

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()