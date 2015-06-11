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

#pragma once

#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::unique_ptr;
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

namespace portcullis { namespace bam {


class BamReader {
        
private:
    
    path bamFile;
    uint16_t threads;
    
    BGZF *fp;
    bam_hdr_t* header;
    bam1_t* c;
    hts_idx_t* index;
    hts_itr_t * iter;
    vector<RefSeq> refs;
    
    BamAlignment b;
    
public:
    
    BamReader(const path& _bamFile, uint16_t threads);
    
    virtual ~BamReader();
    
    vector<RefSeq> getRefs() const { return refs; }
    
    bam_hdr_t* getHeader() const { return header; }
    
    void open();
    
    void close();
    
    bool next();
    
    const BamAlignment& current() const;
        
    void setRegion(const int32_t seqIndex, const int32_t start, const int32_t end);
    
    bool isCoordSortedBam();
};

}}
