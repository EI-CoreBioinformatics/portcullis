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

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <string>
using std::cout;
using std::endl;
using std::stringstream;


#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <portcullis/bam/bam_master.hpp>
#include <portcullis/bam/bam_alignment.hpp>
#include <portcullis/bam/bam_reader.hpp>
#include <portcullis/bam/depth_parser.hpp>
#include <portcullis/bam/genome_mapper.hpp>
using namespace portcullis::bam;

        
TEST(bam, sort) {
    
    string unsortedBam = RESOURCESDIR "/unsorted.bam";
    string sortedBam = RESOURCESDIR "/sorted.test.bam";
    
    string cmd = BamHelper::createSortBamCmd(unsortedBam, sortedBam);
    
    stringstream ss;
    ss << "samtools sort -@ 1 -m 1G " << RESOURCESDIR << "/unsorted.bam " 
            << RESOURCESDIR << "/sorted.test.bam";
    string correct = ss.str();
    
    //cout << "cmd=" << cmd << endl;
    
    // Check the sorted bam file exists
    EXPECT_EQ(cmd, correct);
}


TEST(bam, merge) {
    
    // Merge a couple of BAMs together
    vector<path> bamFiles;
    bamFiles.push_back(RESOURCESDIR "/bam1.bam");
    bamFiles.push_back(RESOURCESDIR "/bam2.bam");
    
    path mergedBam = RESOURCESDIR "/merged.bam";
    string cmd = BamHelper::createMergeBamCmd(bamFiles, mergedBam, 1);
    
    stringstream ss;
    ss << "samtools merge -f -@ 1 " << RESOURCESDIR << "/merged.bam " 
            << RESOURCESDIR << "/bam1.bam " << RESOURCESDIR << "/bam2.bam";
    string correct = ss.str();
    
    EXPECT_EQ(cmd, correct);
}

TEST(bam, is_sorted1) {
    
    string unsortedBam = RESOURCESDIR "/unsorted.bam";
    bool sorted = BamHelper::isCoordSortedBam(unsortedBam);
    
    // Check the merged bam file exists
    EXPECT_EQ(sorted, true);    
}

TEST(bam, is_sorted2) {
    
    string sortedBam = RESOURCESDIR "/sorted.bam";
    bool sorted = BamHelper::isCoordSortedBam(sortedBam);
    
    // Check the merged bam file exists
    EXPECT_EQ(sorted, true);    
}


TEST(bam, depth_test_1) {
    
    DepthParser dp1(RESOURCESDIR "/sorted.bam", 0, true);
    
    vector<uint32_t> batch1;
        
    bool allPos1 = true;
    uint64_t count1 = 0;
    while(dp1.loadNextBatch(batch1)) {

        for(uint32_t cvg : batch1) {
            count1 += cvg;
        }
    }
    
    EXPECT_EQ(allPos1, true);
    //cout << count1 << std::endl;
    
    DepthParser dp2(RESOURCESDIR "/sorted.bam", 0, false);
    
    vector<uint32_t> batch2;
        
    bool allPos2 = true;
    uint64_t count2 = 0;
    while(dp2.loadNextBatch(batch2)) {

        for(uint32_t cvg : batch2) {
            count2 += cvg;
        }        
    }
    
    EXPECT_EQ(allPos2, true);
    //cout << count2 << std::endl;
    
    EXPECT_LE(count2, count1);
}

TEST(bam, genome_mapper_ecoli) {
    
    // Create a new faidx
    bfs::create_directories("temp");
    path in(RESOURCESDIR "/ecoli.fa");
    path out("temp/ecoli.fa");
    
    std::ifstream  src(in.c_str(), std::ios::binary);
    std::ofstream  dst(out.c_str(), std::ios::binary);

    dst << src.rdbuf();
    
    GenomeMapper genomeMapper(out);
    genomeMapper.buildFastaIndex();
    genomeMapper.loadFastaIndex();
    
    // Check faidx file exists
    path faidxFile = genomeMapper.getFastaIndexFile();    
    EXPECT_EQ(bfs::exists(faidxFile), true);
    
    // Check number of seqs is what we expect
    int nbSeqs = genomeMapper.getNbSeqs();    
    EXPECT_EQ(nbSeqs, 1);
    
    // Get seq
    string name = "gi|556503834|ref|NC_000913.3|";
    int len = -1;
    string fullSeq = genomeMapper.fetchBases(name.c_str(), &len);    
    EXPECT_EQ(len, 4641652);
    //BOOST_CHECK(fullSeq != NULL);
    
    string partialSeqExpected = "TCTGACTGCA";
    
    // Get partial seq (method 1 - 1 based)
    string partialSeq1 = genomeMapper.fetchBases((name + ":11-20").c_str(), &len);
    EXPECT_EQ(partialSeqExpected.compare(partialSeq1), 0);
    EXPECT_EQ(len, 10);
    
    // Get partial seq (method 2 - 0 based)
    string partialSeq2 = genomeMapper.fetchBases((char*)name.c_str(), 10, 19, &len);
    EXPECT_EQ(partialSeqExpected.compare(partialSeq2), 0);
    EXPECT_EQ(len, 10);
    
    // Delete the faidx file
    bfs::remove(faidxFile);
}

TEST(bam, padding) {
    
    vector<CigarOp> cigar = CigarOp::createFullCigarFromString("2S14M2I1M1737N8M14S");
    
    string query = "AGAAAGTGGAGAAAAGAATTTGGTGTGGATGATCTTATCACAACCATTCTTTCTGGTGAGACAGAAGC";
    string genomic = "AAAGTGGAGAAAAGAATTTGGTGTGGATGATCTTATCACAACCATTCTTTCTGGTGAGACAGAAGC";
    
    BamAlignment ba;
    ba.setCigar(cigar);
    ba.setRefId(2);
    ba.setPosition(609263);
    ba.setAlignedLength(1787);
    
    int32_t left = 609263;
    int32_t right = 609304;
    string paddedQueryInRegion = ba.getPaddedQuerySeq(query, 609263, 609304, left, right, false);
    string paddedGenomicInRegion = ba.getPaddedGenomeSeq(genomic, 609263, 609304, left, right, false);
    
    EXPECT_EQ(paddedQueryInRegion.size(), paddedGenomicInRegion.size());
    EXPECT_EQ(paddedQueryInRegion, "AAAGTGGAGAAAAGAAT");
    EXPECT_EQ(paddedGenomicInRegion, "AAAGTGGAGAAAAGXXA");
}

TEST(bam, padding2) {
    
    vector<CigarOp> cigar = CigarOp::createFullCigarFromString("14S13M1I2601N9M4918N13M18S");
    
    string query = "ATTGGGGTGTAGATAATTTTATAAAAATTTTTATTTAGGAGGAAAAAAAGGCCGTTTCCAAATATTAC";
    string genomic = "AATTTTATAAAAAAACGGAACTCCGGC";
    
    BamAlignment ba;
    ba.setCigar(cigar);
    ba.setRefId(2);
    ba.setPosition(750577);
    ba.setAlignedLength(7586);
    
    int32_t left = 750577;
    int32_t right = 750603;
    const string paddedQueryInRegion = ba.getPaddedQuerySeq(query, 750577, 750603, left, right, false);
    const string paddedGenomicInRegion = ba.getPaddedGenomeSeq(genomic, 750577, 750603, left, right, false);
    
    EXPECT_EQ(paddedQueryInRegion.size(), paddedGenomicInRegion.size());
    EXPECT_EQ(paddedQueryInRegion, "AATTTTATAAAAAT");
    EXPECT_EQ(paddedGenomicInRegion, "AATTTTATAAAAAX");
}

TEST(bam, padding3) {
    
    vector<CigarOp> cigar = CigarOp::createFullCigarFromString("30S8M25N2M5D28M");
    
    const string query = "ACAAAAACAGAAAAAAAAAGAAAAAAAAATACCAAAACCAACGCCTTCACTTAAAGACAAATATTCAA";
    const string genomic = "TACCAAAG";
    
    BamAlignment ba;
    ba.setCigar(cigar);
    ba.setRefId(2);
    ba.setPosition(4776643);
    ba.setAlignedLength(98);
    
    int32_t left = 4776673;
    int32_t right = 4776680;
    const string paddedQueryInRegion = ba.getPaddedQuerySeq(query, 4776673, 4776680, left, right, false);
    const string paddedGenomicInRegion = ba.getPaddedGenomeSeq(genomic, 4776673, 4776680, left, right, false);
    
    EXPECT_EQ(paddedQueryInRegion.size(), paddedGenomicInRegion.size());
    EXPECT_EQ(paddedQueryInRegion, "CAXXX");
    EXPECT_EQ(paddedGenomicInRegion, "CAAAG");
}
