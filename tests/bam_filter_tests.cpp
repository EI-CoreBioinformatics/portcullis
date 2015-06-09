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
#include <vector>
using std::cout;
using std::endl;
using std::vector;

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "../src/samtools_helper.hpp"
using portcullis::CigarOp;
using portcullis::BamAlignment;

#include "../src/bam_filter.hpp"
using portcullis::BamFilter;


BOOST_AUTO_TEST_SUITE(bam_filter)

BOOST_AUTO_TEST_CASE(completePass) {
    
    vector<CigarOp> completeOkCigar;
    completeOkCigar.push_back(CigarOp('M',10));
    completeOkCigar.push_back(CigarOp('I',2));
    completeOkCigar.push_back(CigarOp('M',10));
    completeOkCigar.push_back(CigarOp('D',10));
    completeOkCigar.push_back(CigarOp('M',10));
    completeOkCigar.push_back(CigarOp('N',10));
    completeOkCigar.push_back(CigarOp('M',10));    
    
    BamAlignment completeOk;
    completeOk.setCigar(completeOkCigar);
    
    //BamFilter filter();
    
    // Check the merged bam file exists
    //BOOST_CHECK(boost::filesystem::exists(mergedBam));    
}

BOOST_AUTO_TEST_SUITE_END()
