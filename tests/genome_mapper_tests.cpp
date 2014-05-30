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

#include <genome_mapper.hpp>

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(genome_mapper)

BOOST_AUTO_TEST_CASE(ecoli)
{
    // Create a new faidx
    portculis::GenomeMapper genomeMapper("resources/ecoli.fa", true, false);
    genomeMapper.buildIndex();
    genomeMapper.loadIndex();
    
    // Check faidx file exists
    string faidxFile = genomeMapper.getFaidxPath();    
    BOOST_CHECK(boost::filesystem::exists(faidxFile));
    
    // Check number of seqs is what we expect
    int nbSeqs = genomeMapper.getNbSeqs();    
    BOOST_CHECK(nbSeqs == 1);
    
    // Get seq
    string name = "gi|556503834|ref|NC_000913.3|";
    int len = -1;
    char* fullSeq = genomeMapper.fetch(name.c_str(), &len);    
    BOOST_CHECK(len == 4641652);
    BOOST_CHECK(fullSeq != NULL);
    
    string partialSeqExpected = "TCTGACTGCA";
    
    // Get partial seq (method 1 - 1 based)
    char* partialSeq1 = genomeMapper.fetch((name + ":11-20").c_str(), &len);
    BOOST_CHECK(partialSeqExpected.compare(partialSeq1) == 0);
    BOOST_CHECK(len == 10);
    
    // Get partial seq (method 2 - 0 based)
    char* partialSeq2 = genomeMapper.fetch((char*)name.c_str(), 10, 19, &len);
    BOOST_CHECK(partialSeqExpected.compare(partialSeq2) == 0);
    BOOST_CHECK(len == 10);
    
    // Delete the faidx file
    boost::filesystem::remove(faidxFile);
}


BOOST_AUTO_TEST_SUITE_END()
