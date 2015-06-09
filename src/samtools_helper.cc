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

#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;
using std::stringstream;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::lexical_cast;

#include <bam.h>

#include "samtools_helper.hpp"
using portcullis::BamAlignment;
using portcullis::BamAlignmentPtr;


// ****** BamReader methods *********

portcullis::BamReader::BamReader(const path& _bamFile, const uint16_t _threads) {
    bamFile = _bamFile;
    threads = _threads;
    
    header = nullptr;    
    index = nullptr;
    iter = nullptr;
    
    c = nullptr;    
}

portcullis::BamReader::~BamReader() {
    
    if (header != nullptr) {
        bam_hdr_destroy(header);
    }
    
    if (c != nullptr) {
        bam_destroy1(c);
    }
        
    if (index != nullptr) {
        hts_idx_destroy(index);
    }
    
    if (iter != nullptr) {
        free(iter->off); 
        free(iter->bins.a); 
        free(iter);
    }
}

void portcullis::BamReader::open() {
    
    // split
    fp = bam_open(bamFile.c_str(), "r");
    if (fp == NULL) {
        BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Could not open input BAM files: ") + bamFile.string()));
    }
    
    // Load header
    header = bam_hdr_read(fp);    
    
    // Identify all the reference sequences in the BAM
    for(uint32_t i = 0; i < header->n_targets; i++) {
        refs.push_back(RefSeq(i, string(header->target_name[i]), header->target_len[i]));
    }
    
    // Load the index
    index = bam_index_load(bamFile.c_str());
    
    // Initialise an empty bam alignment
    c = bam_init1();
}

void portcullis::BamReader::close() {
    bam_close(fp);
}

bool portcullis::BamReader::next() {    
    return bam_iter_read(fp, iter, c) >= 0;    
}

BamAlignmentPtr portcullis::BamReader::current() {
    return make_shared<BamAlignment>(c, false);
}

void portcullis::BamReader::setRegion(const int32_t seqIndex, const int32_t start, const int32_t end) {
    iter = sam_itr_queryi(index, seqIndex, start, end);
}
    

bool portcullis::BamReader::isCoordSortedBam() {    
    string headerText = header->text;   
    return headerText.find("SO:coordinate") != std::string::npos;
}




// ****** Samtools Helper methods *********

path portcullis::SamtoolsHelper::samtoolsExe = "samtools";

bool portcullis::SamtoolsHelper::isCoordSortedBam(const path& bamFile) {

    BGZF *fp;
    
    // split
    fp = strcmp(bamFile.c_str(), "-")? bgzf_open(bamFile.c_str(), "r") : bgzf_dopen(fileno(stdin), "r");
    if (fp == NULL) {
        BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Could not open input BAM files: ") + bamFile.string()));
    }
    
    bam_hdr_t* header = bam_hdr_read(fp);
    
    string headerText = header->text;
   
    bool found = false;
    if (headerText.find("SO:coordinate") != std::string::npos) {
        found = true;        
    }
    
    bam_hdr_destroy(header);
        
    return found;
}

/**
 * Creates a command that can be used to merge multiple BAM files with samtools
 * @param samtoolsExe The path to samtools
 * @param bamFiles The paths to each BAM file to merge
 * @param mergedBamFile The output file
 * @param threads Number of threads to use during merging
 * @return Command line
 */
string portcullis::SamtoolsHelper::createMergeBamCmd(const vector<path>& bamFiles, 
                                                     const path& mergedBamFile, 
                                                     uint16_t threads) {

    stringstream inputFiles;
    for(path p : bamFiles) {            
        inputFiles << " " << p.c_str();
    }

    return samtoolsExe.string() + " merge -f -@ " + lexical_cast<string>(threads) + 
            " " + mergedBamFile.string() + 
            inputFiles.str();
}
    
    
/**
 * Creates a samtools command that can be used to sort a bam file
 * @param samtoolsExe The path to samtools
 * @param unsortedFile The bam file that needs sorting
 * @param sortedFile The path to the new sorted bam file which will be created
 * @param sortByName If true, bam entries are sorted by name, otherwise by position
 * @param threads Number of threads to use
 * @param memory Amount of memory to request
 * @return The command that can be used to sort the bam file
 */
string portcullis::SamtoolsHelper::createSortBamCmd(const path& unsortedFile, 
                                                    const path& sortedFile, 
                                                    bool sortByName, 
                                                    uint16_t threads, 
                                                    const string& memory) {

    return samtoolsExe.string() + " sort -@ " + lexical_cast<string>(threads) + 
            " -m " + memory + " " + (sortByName ? "-n " : "") + unsortedFile.string() + 
            " " + sortedFile.string();
}
    
/**
 * Creates a samtools command that can be used to index a sorted bam file
 * @param samtoolsExe The path to samtools
 * @param sortedBam Path to a sorted bam file to index
 * @return The command that can be used to index the sorted bam file
 */
string portcullis::SamtoolsHelper::createIndexBamCmd(const path& sortedBam) {

    return samtoolsExe.string() + " index " + sortedBam.string();        
}



    