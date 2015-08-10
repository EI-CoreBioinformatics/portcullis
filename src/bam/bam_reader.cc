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

#include <htslib/faidx.h>
#include <bam.h>

#include "bam_alignment.hpp"
using portcullis::bam::BamAlignment;
using portcullis::bam::BamAlignmentPtr;
using portcullis::bam::RefSeq;

#include "bam_reader.hpp"

// ****** BamReader methods *********

portcullis::bam::BamReader::BamReader(const path& _bamFile, const uint16_t _threads) {
    bamFile = _bamFile;
    threads = _threads;
    
    header = nullptr;    
    index = nullptr;
    iter = nullptr;
    
    c = nullptr;    
}

portcullis::bam::BamReader::~BamReader() {
    
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

void portcullis::bam::BamReader::open() {
    
    // split
    fp = bam_open(bamFile.c_str(), "r");
    if (fp == NULL) {
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
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
    b.setRaw(c);
}

void portcullis::bam::BamReader::close() {
    bam_close(fp);
}

const unordered_map<int32_t, RefSeq>& portcullis::bam::BamReader::calcRefMap() {
    
    refMap.clear();
    
    for(auto& ref : refs) {
        refMap[ref.index] = ref;
    }
    
    return refMap;
}

bool portcullis::bam::BamReader::next() {
    bool res = bam_iter_read(fp, iter, c) >= 0;
    b.setRaw(c);    
    return res;    
}

const BamAlignment& portcullis::bam::BamReader::current() const {
    return b;
}

void portcullis::bam::BamReader::setRegion(const int32_t seqIndex, const int32_t start, const int32_t end) {
    iter = sam_itr_queryi(index, seqIndex, start, end);
}
    

bool portcullis::bam::BamReader::isCoordSortedBam() {    
    string headerText = header->text;   
    return headerText.find("SO:coordinate") != std::string::npos;
}

