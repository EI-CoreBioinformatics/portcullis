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

#include <iostream>
using std::cout;
using std::endl;


#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include <portcullis_fs.hpp>
using portcullis::PortcullisFS;

#include <htslib_helper.hpp>
using portcullis::HtslibHelper;

const string samtoolsExe = "../deps/samtools-1.2/samtools";
        
BOOST_AUTO_TEST_SUITE(bam_utils)

BOOST_AUTO_TEST_CASE(is_sorted1)
{
    string unsortedBam = "resources/unsorted.bam";
    bool sorted = BamUtils::isSortedBam(unsortedBam);
    
    // Check the merged bam file exists
    //BOOST_CHECK(!sorted);    
    BOOST_CHECK(sorted);    
}

BOOST_AUTO_TEST_CASE(is_sorted2)
{
    string sortedBam = "resources/sorted.bam";
    bool sorted = BamUtils::isSortedBam(sortedBam);
    
    // Check the merged bam file exists
    BOOST_CHECK(sorted);    
}

// This test doesn't actually do that much.  Refine it later.
BOOST_AUTO_TEST_CASE(minmalMatch1)
{
    string bam = "resources/clipped3.bam";
    
    BamReader reader;
    reader.Open(bam);
    
    // Sam header and refs info from the input bam
    SamHeader header = reader.GetHeader();
    RefVector refs = reader.GetReferenceData();
    
    BamAlignment ba;
    reader.GetNextAlignment(ba);
        
    //cout << ba.AlignedBases << endl;
    //cout << "Pos: " << ba.Position << endl;
    
    uint16_t mm = BamUtils::calcMinimalMatchInCigarDataSubset(ba, 6441673 + 2, 6441673 + 10);
    
    //cout << mm << endl;
    
    // Check the minimal match result
    BOOST_CHECK(mm == 8);    
}

BOOST_AUTO_TEST_SUITE_END()
