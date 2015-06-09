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
#include <boost/algorithm/string.hpp>

#include <samtools_helper.hpp>
using portcullis::SamtoolsHelper;
using portcullis::BamReader;
using portcullis::BamAlignment;
using portcullis::BamAlignmentPtr;


        
BOOST_AUTO_TEST_SUITE(samtools_helper)

BOOST_AUTO_TEST_CASE(sort) {
    
    SamtoolsHelper::samtoolsExe = "../deps/samtools-1.2/samtools";

    string unsortedBam = "resources/unsorted.bam";
    string sortedBam = "resources/sorted.test.bam";
    
    string cmd = SamtoolsHelper::createSortBamCmd(unsortedBam, sortedBam);
    
    string correct("../deps/samtools-1.2/samtools sort -@ 1 -m 1G resources/unsorted.bam resources/sorted.test.bam");
    
    //cout << "cmd=" << cmd << endl;
    
    // Check the sorted bam file exists
    BOOST_CHECK(boost::equals(cmd, correct));
}


BOOST_AUTO_TEST_CASE(merge) {
    
    SamtoolsHelper::samtoolsExe = "../deps/samtools-1.2/samtools";

    // Merge a couple of BAMs together
    vector<path> bamFiles;
    bamFiles.push_back("resources/bam1.bam");
    bamFiles.push_back("resources/bam2.bam");
    
    path mergedBam = "resources/merged.bam";
    string cmd = SamtoolsHelper::createMergeBamCmd(bamFiles, mergedBam, 1);
    
    string correct("../deps/samtools-1.2/samtools merge -f -@ 1 resources/merged.bam resources/bam1.bam resources/bam2.bam");
    
    BOOST_CHECK(boost::equals(cmd, correct));
}

BOOST_AUTO_TEST_CASE(is_sorted1)
{
    string unsortedBam = "resources/unsorted.bam";
    bool sorted = SamtoolsHelper::isCoordSortedBam(unsortedBam);
    
    // Check the merged bam file exists
    //BOOST_CHECK(!sorted);    
    BOOST_CHECK(sorted);    
}

BOOST_AUTO_TEST_CASE(is_sorted2)
{
    string sortedBam = "resources/sorted.bam";
    bool sorted = SamtoolsHelper::isCoordSortedBam(sortedBam);
    
    // Check the merged bam file exists
    BOOST_CHECK(sorted);    
}

// This test doesn't actually do that much.  Refine it later.
BOOST_AUTO_TEST_CASE(minmalMatch1)
{
    path bam = "resources/clipped3.bam";
    
    BamReader reader(bam, 1);
    reader.open();
    
    // Sam header and refs info from the input bam
    reader.next();
    BamAlignmentPtr ba = reader.current();
        
    //cout << ba.AlignedBases << endl;
    //cout << "Pos: " << ba.Position << endl;
    
    uint16_t mm = ba->calcMinimalMatchInCigarDataSubset(6441673 + 2, 6441673 + 10);
    
    //cout << mm << endl;
    
    // Check the minimal match result
    BOOST_CHECK(mm == 8);    
}

BOOST_AUTO_TEST_CASE(depth_test_1)
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

BOOST_AUTO_TEST_CASE(genome_mapper_ecoli)
{
    // Create a new faidx
    portcullis::GenomeMapper genomeMapper("resources/ecoli.fa");
    genomeMapper.buildFastaIndex();
    genomeMapper.loadFastaIndex();
    
    // Check faidx file exists
    path faidxFile = genomeMapper.getFastaIndexFile();    
    BOOST_CHECK(boost::filesystem::exists(faidxFile));
    
    // Check number of seqs is what we expect
    int nbSeqs = genomeMapper.getNbSeqs();    
    BOOST_CHECK(nbSeqs == 1);
    
    // Get seq
    string name = "gi|556503834|ref|NC_000913.3|";
    int len = -1;
    char* fullSeq = genomeMapper.fetchBases(name.c_str(), &len);    
    BOOST_CHECK(len == 4641652);
    BOOST_CHECK(fullSeq != NULL);
    
    string partialSeqExpected = "TCTGACTGCA";
    
    // Get partial seq (method 1 - 1 based)
    char* partialSeq1 = genomeMapper.fetchBases((name + ":11-20").c_str(), &len);
    BOOST_CHECK(partialSeqExpected.compare(partialSeq1) == 0);
    BOOST_CHECK(len == 10);
    
    // Get partial seq (method 2 - 0 based)
    char* partialSeq2 = genomeMapper.fetchBases((char*)name.c_str(), 10, 19, &len);
    BOOST_CHECK(partialSeqExpected.compare(partialSeq2) == 0);
    BOOST_CHECK(len == 10);
    
    // Delete the faidx file
    boost::filesystem::remove(faidxFile);
}

BOOST_AUTO_TEST_SUITE_END()