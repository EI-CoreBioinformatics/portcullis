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
using std::shared_ptr;
using std::make_shared;
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
#include <htslib/sam.h>
#include <htslib/bgzf.h>

namespace portcullis { namespace bam {


typedef struct {     // auxiliary data structure
	BGZF* fp;      // the file handler
	hts_itr_t* iter; // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

typedef struct {
    int ref;
    int32_t pos;
    uint32_t depth;    
} depth;

class DepthParser {
private:
 
    // Path to the original genome file in fasta format
    path bamFile;
    uint8_t strandSpecific;
    bool allowGappedAlignments;
    
    bam_hdr_t *header;
    aux_t** data;
    bam_mplp_t mplp;
        
    depth last;
    bool start;
    int res;
    
protected:
    
    // This function reads a BAM alignment from one BAM file.
    static int read_bam(void *data, bam1_t *b);
    
    // This function reads a BAM alignment from one BAM file.
    static int read_bam_skip_gapped(void *data, bam1_t *b);
    
    
public:
    
    DepthParser(path _bamFile, uint8_t _strandSpecific, bool _allowGappedAlignments);
    
    virtual ~DepthParser();
    
     
    string getCurrentRefName() const {
        return string(header->target_name[last.ref]);
    }
    
    int32_t getCurrentRefIndex() const {
        return last.ref;
    }
    
    

    bool loadNextBatch(vector<uint32_t>& depths);
    
};

}}
