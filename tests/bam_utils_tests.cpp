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

#include <bam_utils.hpp>
#include <bamtools_sort.cpp>

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(bam_utils)

BOOST_AUTO_TEST_CASE(merge)
{
    // Merge a couple of BAMs together
    vector<string> bamFiles;
    bamFiles.push_back("resources/bam1.bam");
    bamFiles.push_back("resources/bam2.bam");
    
    string mergedBam = "resources/merged.bam";
    portculis::mergeBams(bamFiles, mergedBam);
    
    // Check the merged bam file exists
    BOOST_CHECK(boost::filesystem::exists(mergedBam));
    
    // Delete the merged bam file
    boost::filesystem::remove(mergedBam);
}

BOOST_AUTO_TEST_CASE(sort)
{
    string unsortedBam = "resources/unsorted.bam";
    string sortedBam = "resources/sorted.test.bam";
    portculis::sortBam(unsortedBam, sortedBam, false);
    
    // Check the sorted bam file exists
    BOOST_CHECK(boost::filesystem::exists(sortedBam));
    
    // Check the sorted bam file is actually sorted
    BOOST_CHECK(portculis::isSortedBam(sortedBam)); 
    
    // Delete the sorted bam file
    boost::filesystem::remove(sortedBam);
}

BOOST_AUTO_TEST_CASE(is_sorted1)
{
    string unsortedBam = "resources/unsorted.bam";
    bool sorted = portculis::isSortedBam(unsortedBam);
    
    // Check the merged bam file exists
    //BOOST_CHECK(!sorted);    
    BOOST_CHECK(sorted);    
}

BOOST_AUTO_TEST_CASE(is_sorted2)
{
    string sortedBam = "resources/sorted.bam";
    bool sorted = portculis::isSortedBam(sortedBam);
    
    // Check the merged bam file exists
    BOOST_CHECK(sorted);    
}


BOOST_AUTO_TEST_CASE(index)
{
    string sortedBam = "resources/sorted.bam";
    string indexedBam = "resources/sorted.bam.bti";
    portculis::indexBam(sortedBam);
    
    // Check the indexed bam file exists
    BOOST_CHECK(boost::filesystem::exists(indexedBam));
    
    // Delete the indexed bam file
    boost::filesystem::remove(indexedBam);
}

BOOST_AUTO_TEST_SUITE_END()
