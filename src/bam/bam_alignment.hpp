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
#include <bam.h>

#include "bam_master.hpp"

namespace portcullis { namespace bam {

// CIGAR constants
const char BAM_CIGAR_MATCH_CHAR    = 'M';
const char BAM_CIGAR_INS_CHAR      = 'I';
const char BAM_CIGAR_DEL_CHAR      = 'D';
const char BAM_CIGAR_REFSKIP_CHAR  = 'N';
const char BAM_CIGAR_SOFTCLIP_CHAR = 'S';
const char BAM_CIGAR_HARDCLIP_CHAR = 'H';
const char BAM_CIGAR_PAD_CHAR      = 'P';
const char BAM_CIGAR_EQUAL_CHAR    = '=';
const char BAM_CIGAR_DIFF_CHAR     = 'X';
const char BAM_CIGAR_BACK_CHAR     = 'B';


struct CigarOp {
    char     type;
    uint32_t length;
    
    CigarOp(char _type, uint32_t _length) {
        type = _type;
        length = _length;
    }
    
    CigarOp(const string& cigar);
    
    string toString() const {
        return lexical_cast<string>(length) + type;
    }
    
    static vector<CigarOp> createFullCigarFromString(const string& cigar);
    
    static inline bool opConsumesQuery(char op) {
        
        switch (op) {
            case BAM_CIGAR_MATCH_CHAR :
            case BAM_CIGAR_INS_CHAR :
            case BAM_CIGAR_SOFTCLIP_CHAR :
            case BAM_CIGAR_EQUAL_CHAR :
            case BAM_CIGAR_DIFF_CHAR :
                return true;
            default:
                return false;
        }
    }
    
    static inline bool opConsumesReference(char op) {
        
        switch (op) {
            case BAM_CIGAR_MATCH_CHAR :
            case BAM_CIGAR_DEL_CHAR :            
            case BAM_CIGAR_REFSKIP_CHAR :
            case BAM_CIGAR_EQUAL_CHAR :
            case BAM_CIGAR_DIFF_CHAR :
                return true;
            default:
                return false;
        }
    }
};

class BamAlignment {
    
private:
    
    bam1_t* b;
    bool managed;
    
    uint32_t alFlag;
    int32_t position;
    int32_t alignedLength;
    int32_t refId;
    vector<CigarOp> cigar;
    
    void init();
        
public:

    /**
     * Creates an empty samtools bam alignment
     */    
    BamAlignment();
    
    /**
     * Makes a deep copy of the provided samtools bam alignment
     * @param _b Samtools bam alignment
     */
    BamAlignment(bam1_t* _b, bool duplicate);
    
    /**
     * Makes a deep copy of an existing BamAlignment
     * @param other
     */
    BamAlignment(const shared_ptr<BamAlignment> other);
    
    /**
     * Makes a deep copy of an existing BamAlignment
     * @param other
     */
    BamAlignment(const BamAlignment& other);
    
    /**
     * Deletes the underlying samtools bam alignment only if it is managed (owned) 
     * by this object
     */
    virtual ~BamAlignment();
    
    void setRaw(bam1_t* b);
    
    bam1_t* getRaw() const;
     
    void setCigar(vector<CigarOp>& cig) {
        cigar = cig;
    }
    
    const vector<CigarOp>& getCigar() const {
        return cigar;
    }
    
    const string getCigarAsString() const {
        
        stringstream ss;
        for(const auto& c : cigar) {
            ss << c.toString();
        }
        return ss.str();
    }
    
    void setCigarOpAt(uint32_t index, CigarOp cigarOp) {
        cigar[index] = cigarOp;
    }
    
    void setAlignedLength(int32_t alignedLength) {
        this->alignedLength = alignedLength;
    }

    void setPosition(int32_t position) {
        this->position = position;
    }

    void setRefId(int32_t refId) {
        this->refId = refId;
    }

    
    const CigarOp& getCigarOpAt(uint32_t index) const {
        return cigar[index];
    }
    
    size_t getNbCigarOps() const {
        return cigar.size();
    }
    
    int32_t getPosition() const {
        return position;
    }
    
    int32_t getStart() const {
        return position;
    }
    
    int32_t getStart(bool afterClipping) const {
        return afterClipping && cigar.front().type == BAM_CIGAR_SOFTCLIP_CHAR ? position + cigar.front().length : position;    
    }
    
    int32_t getEnd() const {
        return position + alignedLength;
    }
    
    int32_t getEnd(bool afterClipping) const {
        return afterClipping && cigar.back().type == BAM_CIGAR_SOFTCLIP_CHAR ? getEnd() - cigar.back().length : getEnd();    
    }
    
    int32_t getReferenceId() const {
        return refId;
    }
    
    int32_t getLength() const {
        return b->core.l_qseq;
    }
    
    int32_t getMapQuality() const {
        return b->core.qual;
    }
    
    bool isDuplicate() const {
        return (alFlag & BAM_FDUP) != 0;
    }

    bool isFailedQC() const {
        return (alFlag & BAM_FQCFAIL) != 0;
    }

    bool isFirstMate() const {
        return (alFlag & BAM_FREAD1) != 0;
    }

    bool isMapped() const {
        return (alFlag & BAM_FUNMAP) == 0;
    }

    bool isMateMapped() const {
        return (alFlag & BAM_FMUNMAP) == 0;
    }

    bool isMateReverseStrand() const {
        return (alFlag & BAM_FMREVERSE) != 0;
    }

    bool isPaired() const {
        return (alFlag & BAM_FPAIRED) != 0;
    }

    bool isPrimaryAlignment() const  {
        return (alFlag & BAM_FSECONDARY) == 0;
    }

    bool isProperPair() const {
        return (alFlag & BAM_FPROPER_PAIR) != 0;
    }

    bool isReverseStrand() const {
        return (alFlag & BAM_FREVERSE) != 0;
    }

    bool isSecondMate() const {
        return (alFlag & BAM_FREAD2) != 0;
    }
    
    
    
    string deriveName() const;
    
    string getQuerySeq() const;
    string getQuerySeqAfterClipping() const;
    string getQuerySeqAfterClipping(const string& query_seq) const;
    
    bool isSplicedRead() const;
    
    uint32_t getNbJunctionsInRead() const;
    
    bool isMultiplySplicedRead() const {        
        return getNbJunctionsInRead() > 1;
    }
    
    uint32_t calcNbAlignedBases(int32_t start, int32_t end, bool includeSoftClips) const;
    
    string getPaddedQuerySeq(uint32_t start, uint32_t end, uint32_t& actual_start, uint32_t& actual_end, const bool include_soft_clips) const;
    string getPaddedQuerySeq(const string& querySeq, uint32_t start, uint32_t end, uint32_t& actual_start, uint32_t& actual_end, const bool include_soft_clips) const;
    string getPaddedGenomeSeq(const string& fullGenomeSeq, uint32_t start, uint32_t end, uint32_t q_start, uint32_t q_end, const bool include_soft_clips) const;
    
    string toString() const;
    string toString(bool afterClipping) const;
    
};

typedef shared_ptr<BamAlignment> BamAlignmentPtr;

}}
