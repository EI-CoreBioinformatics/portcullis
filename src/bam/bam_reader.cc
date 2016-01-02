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
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using std::endl;
using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;
using std::stringstream;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>  
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::lexical_cast;

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "bam_master.hpp"
#include "bam_alignment.hpp"
using portcullis::bam::BamAlignment;
using portcullis::bam::BamAlignmentPtr;

#include "bam_reader.hpp"

// ****** BamReader methods *********

portcullis::bam::BamReader::BamReader(const path& _bamFile) {
    bamFile = _bamFile;
    
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
    fp = bgzf_open(bamFile.c_str(), "r");
    if (fp == NULL) {
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Could not open input BAM files: ") + bamFile.string()));
    }
    
    // Load header
    header = bam_hdr_read(fp);
    
    // Load the index
    index = bam_index_load(bamFile.c_str());
    
    // Initialise an empty bam alignment
    c = bam_init1();
    b.setRaw(c);
}

void portcullis::bam::BamReader::close() {
    bgzf_close(fp);
}


/**
 * Creates a list of the reference target sequences stored in this BAM
 * @return 
 */
shared_ptr<RefSeqPtrList> portcullis::bam::BamReader::createRefList() {
    
    shared_ptr<RefSeqPtrList> refs = make_shared<RefSeqPtrList>();
    
    // Identify all the reference sequences in the BAM
    for(uint32_t i = 0; i < header->n_targets; i++) {
        refs->push_back(make_shared<RefSeq>(i, string(header->target_name[i]), header->target_len[i]));
    }
    
    return refs;    
}


/**
 * Makes a map directly from the BamReader, keyed off the reference index
 * @return 
 */
shared_ptr<RefSeqPtrIndexMap> portcullis::bam::BamReader::createRefMap() {    
    return createRefMap(*createRefList());
}

/**
 * Makes a map based on the same RefSeqs stored in the provided list, keyed off the reference index
 * @param list
 * @return 
 */
shared_ptr<RefSeqPtrIndexMap> portcullis::bam::BamReader::createRefMap(const RefSeqPtrList& list) {
    
    shared_ptr<RefSeqPtrIndexMap> refMap = make_shared<RefSeqPtrIndexMap>();
    
    for(auto& ref : list) {
        (*refMap)[ref->index] = ref;
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

string portcullis::bam::BamReader::bamDetails() const {
    
    stringstream ss;

    string headerText = this->header->text;
    
    vector<string> lines;
    boost::split( lines, headerText, boost::is_any_of("\n"), boost::token_compress_on );
    
    ss << "BAM details:" << endl
         << " - File path: " << bamFile << endl
         << " - # Target sequences: " << this->header->n_targets << endl;
         
    string hd = "@HD";
    string pg = "@PG";
    uint16_t programCount = 1;
    for(auto& line : lines) {
        if (line.compare(0, hd.length(), hd) == 0) {
            
            string version = "VN:";
            string sortOrder = "SO:";
            string ag = "GO:";
            vector<string> parts;
            boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );
            for (auto& part : parts) {
                if (part.compare(0, version.length(), version) == 0) {
                    ss << " - Format version: " << part.substr(version.length()) << endl;
                }
                else if (part.compare(0, sortOrder.length(), sortOrder) == 0) {
                    ss << " - Sort order: " << part.substr(sortOrder.length()) << endl;
                }
                else if (part.compare(0, ag.length(), ag) == 0) {
                    ss << " - Alignment grouping: " << part.substr(ag.length()) << endl;
                }
            }            
        }
        else if (line.compare(0, pg.length(), pg) == 0) {
            ss << " - Program " << programCount++ << ": " << endl;
            string id = "ID:";
            string name = "PN:";
            string cl = "CL:";
            string vn = "VN:";
            string desc = "DS:";
            vector<string> parts;
            boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );
            for (auto& part : parts) {
                if (part.compare(0, id.length(), id) == 0) {
                    ss << "  - ID: " << part.substr(id.length()) << endl;
                }
                else if (part.compare(0, name.length(), name) == 0) {
                    ss << "  - Name: " << part.substr(name.length()) << endl;
                }
                else if (part.compare(0, cl.length(), cl) == 0) {
                    ss << "  - Command line: " << part.substr(cl.length()) << endl;
                }
                else if (part.compare(0, vn.length(), vn) == 0) {
                    ss << "  - Version: " << part.substr(vn.length()) << endl;
                }
                else if (part.compare(0, desc.length(), desc) == 0) {
                    ss << "  - Description: " << part.substr(desc.length()) << endl;
                }
            }
        }
    }
    
    return ss.str();
}
        

