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

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::cout;
using std::endl;
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
#include <bam.h>

#include "bam_alignment.hpp"
    
void portcullis::bam::BamAlignment::init() {
    alFlag = b->core.flag;
    position = b->core.pos;
    refId = b->core.tid;
    uint32_t* c = bam_get_cigar(b);
    alignedLength = 0;

    cigar.clear();
    for(uint32_t i = 0; i < b->core.n_cigar; i++) {
        CigarOp op(bam_cigar_opchr(c[i]), bam_cigar_oplen(c[i]));
        cigar.push_back(op);
        if (CigarOp::opConsumesReference(op.type)) {
            alignedLength += op.length;
        }
    }
    
}
   
/**
 * Creates an empty samtools bam alignment
 */    
portcullis::bam::BamAlignment::BamAlignment() {
    b = bam_init1();
    managed = true;
    alFlag = 0;
    position = -1;
    alignedLength = 0;
    refId = -1;
}

/**
 * Makes a deep copy of the provided samtools bam alignment
 * @param _b Samtools bam alignment
 */
portcullis::bam::BamAlignment::BamAlignment(bam1_t* _b, bool duplicate) {
    b = duplicate ? bam_dup1(_b) : _b;
    managed = duplicate;
    init();
}

/**
 * Makes a deep copy of an existing BamAlignment
 * @param other
 */
portcullis::bam::BamAlignment::BamAlignment(const shared_ptr<BamAlignment> other) {
    b = bam_dup1(other->b);
    managed = true;
    init();
}

/**
 * Makes a deep copy of an existing BamAlignment
 * @param other
 */
portcullis::bam::BamAlignment::BamAlignment(const BamAlignment& other) {
    b = bam_dup1(other.b);
    managed = true;
    init();
}

/**
 * Deletes the underlying samtools bam alignment only if it is managed (owned) 
 * by this object
 */
portcullis::bam::BamAlignment::~BamAlignment() {
    if (managed) 
        bam_destroy1(b);
}

void portcullis::bam::BamAlignment::setRaw(bam1_t* b) {
    this->b = b;
    managed = false;
    init();
}

bam1_t* portcullis::bam::BamAlignment::getRaw() const {
    return b;
}



string portcullis::bam::BamAlignment::deriveName() const {

    string qName = bam_get_qname(b);

    return isPaired() ?
            qName + (this->isFirstMate() ? 
                        "_R1" :
                        this->isSecondMate() ? 
                            "_R2" :
                            "_R?") :
            qName;
}

string portcullis::bam::BamAlignment::getQuerySeq() const {
    
    stringstream ss;
    for (uint32_t i = 0; i < b->core.l_qseq; ++i) {
       ss << seq_nt16_str[bam_seqi(bam_get_seq(b),i)];
    }
    
    return ss.str();
}

bool portcullis::bam::BamAlignment::isSplicedRead() const {
    for(const auto& op : cigar) {
        if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
            return true;
        }
    }

    return false;
}

uint32_t portcullis::bam::BamAlignment::getNbJunctionsInRead() const {

    int32_t nbJunctions = 0;
    for(const auto& op : cigar) {
        if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
            nbJunctions++;
        }
    }

    return nbJunctions;
}

uint32_t portcullis::bam::BamAlignment::calcNbAlignedBases(int32_t start, int32_t end) const {

    if (start > getEnd() || end < position) {
        string align = this->toString();
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Found an alignment that does not have a presence in the requested region.  Requested region: ") +
                lexical_cast<string>(start) + "-" + lexical_cast<string>(end) + ".  Alignment: " + align));
    }

    int32_t count = 0;
    int32_t pos = position;
    for(const auto& op : cigar) {

        if (pos > end) {
            break;
        }
        
        if (CigarOp::opConsumesReference(op.type)) {
            if (pos >= start) {
                count += op.length;
            }
            pos += op.length;
        }
    } 

    return count;
}


string portcullis::bam::BamAlignment::getPaddedQuerySeq(uint32_t start, uint32_t end, uint32_t& actual_start, uint32_t& actual_end) const {
    
    if (start > getEnd() || end < position)
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Found an alignment that does not have a presence in the requested region")));

    int32_t qPos = 0;
    int32_t rPos = position;
    
    string query = this->getQuerySeq();
    stringstream ss;
    
    for(const auto& op : cigar) {

        bool consumesRef = CigarOp::opConsumesReference(op.type);
        bool consumesQuery = CigarOp::opConsumesQuery(op.type);
        
        // Skips any cigar ops before start position
        if (rPos < start) {
            if (consumesRef) rPos += op.length;
            if (consumesQuery) qPos += op.length;            
            continue;
        }
        
        // Stop once we get to the end of the region, and make sure we don't end on a refskip
        if (rPos >= end || (op.type == BAM_CIGAR_REFSKIP_CHAR && rPos + op.length >= end)) break;
        
        if (consumesQuery && op.type != BAM_CIGAR_SOFTCLIP_CHAR) {
            
            // Don't return anything that runs off the end cap
            int32_t len = rPos + op.length > end ? end - rPos + 1 : op.length;
            
            //cout << qPos << endl;
            if (qPos < 0 || qPos + len > query.size()) {
                BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                    "Can't extract cigar op sequence from query string.  Current position in query: ") 
                        + lexical_cast<string>(qPos) + " Cigar op: " + op.type + lexical_cast<string>(op.length)));
            }
            
            ss << query.substr(qPos, len);
        }
        else if (consumesRef) { // i.e. consumes reference but not query (DEL or REF_SKIP ops)
            string s;
            s.resize(op.length);
            std::fill(s.begin(), s.end(), 'N'); 
            ss << s;
        }
        
        if (consumesRef) rPos += op.length;
        if (consumesQuery) qPos += op.length;
    }
    
    actual_start = position > start ? position : start;
    actual_end = rPos < end ? rPos : end;
    
    return ss.str();
}


string portcullis::bam::BamAlignment::getPaddedGenomeSeq(const string& genomeSeq, uint32_t start, uint32_t end, uint32_t q_start, uint32_t q_end) const {
    
    if (start > getEnd() || end < position)
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Found an alignment that does not have a presence in the requested region")));

    //int32_t pos = 0;
    //int32_t qPos = 0;
    int32_t rPos = position;
    int32_t startDelta = q_start - start;
    int32_t endDelta = end - q_end;
    
    if (startDelta < 0) {
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Query start position was before genomic region start position.  Query start: ") + lexical_cast<string>(q_start) + "; Genomic start: " + lexical_cast<string>(start)));
    }
    
    if (endDelta < 0) {
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Query end position was beyond genomic region end position.  Query end: ") + lexical_cast<string>(q_end) + "; Genomic end: " + lexical_cast<string>(end)));
    }
    
    
    stringstream ss;
        
    for(const auto& op : cigar) {

        bool consumesRef = CigarOp::opConsumesReference(op.type);
        bool consumesQuery = CigarOp::opConsumesQuery(op.type);
        
        // Skips any cigar ops before start position
        if (rPos < q_start) {
            if (consumesRef) rPos += op.length;
            //if (consumesQuery) qPos += op.length;
            continue;
        }
        
        // Ends cigar loop once we reach the end of the region
        if (rPos >= q_end) break;
        
        if (consumesRef) {
            
            int32_t seqOffset = rPos - start;
            int32_t len = rPos + op.length > q_end ? q_end - rPos + 1 : op.length;
            
            if (seqOffset < 0 || seqOffset + len > genomeSeq.size()) {
                BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                    "Can't extract cigar op sequence from extracted genome region.\nCurrent position in extracted genome region: ") 
                        + lexical_cast<string>(seqOffset) + 
                        "; Genome region length: " + lexical_cast<string>(genomeSeq.size()) +
                        "\nFull cigar: " + this->getCigarAsString() +
                        "; Current cigar op: " + op.type + lexical_cast<string>(op.length) +
                        "\nCurrent genomic position: " + lexical_cast<string>(rPos) +
                        "\nAlignment region: " + lexical_cast<string>(position) + "," + lexical_cast<string>(position + alignedLength) +
                        "\nGenome region: " + lexical_cast<string>(start) + "," + lexical_cast<string>(end) +
                        "\nQuery region: " + lexical_cast<string>(q_start) + "," + lexical_cast<string>(q_end)));
            }
            
            //cout << seqOffset << endl;
            ss << genomeSeq.substr(seqOffset, len);
        }
        else if (consumesQuery && op.type != BAM_CIGAR_SOFTCLIP_CHAR) {   // Consumes query but not reference ('I' op)
            string s;
            s.resize(op.length);
            std::fill(s.begin(), s.end(), 'N');
            ss << s;
        }     
        
        if (consumesRef) rPos += op.length;
        //if (consumesQuery) qPos += op.length;
    } 

    return ss.str();
}

string portcullis::bam::BamAlignment::toString() const {
    
    stringstream ss;
    ss << refId << "(" << position << "-" << getEnd() << ")";
    return ss.str();
}