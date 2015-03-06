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

#include <bam_utils.hpp>

using std::cout;
using std::endl;

using portcullis::bamtools::BamUtils;

BOOST_AUTO_TEST_SUITE(bam_utils)

BOOST_AUTO_TEST_CASE(merge) {
    
    // Merge a couple of BAMs together
    vector<string> bamFiles;
    bamFiles.push_back("resources/bam1.bam");
    bamFiles.push_back("resources/bam2.bam");
    
    string mergedBam = "resources/merged.bam";
    BamUtils::mergeBams(bamFiles, mergedBam);
    
    // Check the merged bam file exists
    BOOST_CHECK(boost::filesystem::exists(mergedBam));
    
    // Delete the merged bam file
    boost::filesystem::remove(mergedBam);
}

BOOST_AUTO_TEST_CASE(sort) {
    
    string unsortedBam = "resources/unsorted.bam";
    string sortedBam = "resources/sorted.test.bam";
    
    string cmd = BamUtils::createSortBamCmd(unsortedBam, sortedBam);
    
    string correct("samtools sort -@ 1 -m 1G resources/unsorted.bam resources/sorted.test.bam");
    
    //cout << "cmd=" << cmd << endl;
    
    // Check the sorted bam file exists
    BOOST_CHECK(cmd == correct);
}

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
    
    // Check the merged bam file exists
    BOOST_CHECK(mm == 0);    
}

BOOST_AUTO_TEST_SUITE_END()
