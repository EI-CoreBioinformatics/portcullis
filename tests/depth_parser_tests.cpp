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

#include <boost/filesystem.hpp>

#include <depth_parser.hpp>
using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(depth_parser)

BOOST_AUTO_TEST_CASE(test1)
{
    portcullis::DepthParser dp1("resources/sorted.bam", 0, true);
    
    vector<uint32_t> batch1;
        
    bool allPos1 = true;
    uint64_t count1 = 0;
    while(dp1.loadNextBatch(batch1)) {

        for(uint32_t cvg : batch1) {
            if (cvg < 0) {
                allPos1 = false;                
            }
            else {
                count1 += cvg;
            }
        }
    }
    
    BOOST_CHECK(allPos1);
    //cout << count1 << std::endl;
    
    portcullis::DepthParser dp2("resources/sorted.bam", 0, false);
    
    vector<uint32_t> batch2;
        
    bool allPos2 = true;
    uint64_t count2 = 0;
    while(dp2.loadNextBatch(batch2)) {

        for(uint32_t cvg : batch2) {
            if (cvg < 0) {
                allPos2 = false;                
            }
            else {
                count2 += cvg;
            }
        }        
    }
    
    BOOST_CHECK(allPos2);
    //cout << count2 << std::endl;
    
    BOOST_CHECK(count2 <= count1);
}


BOOST_AUTO_TEST_SUITE_END()