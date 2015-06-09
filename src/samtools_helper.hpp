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

#include <bam.h>

namespace portcullis {

typedef boost::error_info<struct SamtoolsError,string> SamtoolsErrorInfo;
struct SamtoolsException: virtual boost::exception, virtual std::exception { };

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
    int32_t refId;
    vector<CigarOp> cigar;
    
    void init() {
        alFlag = b->core.flag;
        position = b->core.pos;
        refId = b->core.tid;
        uint32_t* c = bam_get_cigar(b);
        
        for(uint32_t i = 0; i < b->core.n_cigar; i++) {
            cigar.push_back(CigarOp(bam_cigar_opchr(c[i]), bam_cigar_oplen(c[i])));
        }
    }
    
public:

    /**
     * Creates an empty samtools bam alignment
     */    
    BamAlignment() {
        b = bam_init1();
        managed = true;
        alFlag = 0;
        position = -1;
        refId = -1;
    }
    
    /**
     * Makes a deep copy of the provided samtools bam alignment
     * @param _b Samtools bam alignment
     */
    BamAlignment(bam1_t* _b, bool duplicate) {
        b = duplicate ? bam_dup1(_b) : _b;
        managed = duplicate;
        init();
    }
    
    /**
     * Makes a deep copy of an existing BamAlignment
     * @param other
     */
    BamAlignment(const BamAlignment& other) {
        b = bam_dup1(other.b);
        managed = true;
        init();
    }
    
    /**
     * Deletes the underlying samtools bam alignment only if it is managed (owned) 
     * by this object
     */
    virtual ~BamAlignment() {
        if (managed) bam_destroy1(b);
    }
    
    bam1_t* getRaw() const {
        return b;
    }
        
    vector<CigarOp> getCigar() const {
        return cigar;
    }
    
    void setCigarOpAt(uint32_t index, CigarOp cigarOp) {
        cigar[index] = cigarOp;
    }
    
    CigarOp getCigarOpAt(uint32_t index) const {
        return cigar[index];
    }
    
    size_t getNbCigarOps() const {
        return cigar.size();
    }
    
    int32_t getPosition() const {
        return position;
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
    
    string deriveName() const {
        
        string qName = bam_get_qname(b);
        
        return this->isPaired() ?
                qName + (this->isFirstMate() ? 
                            "_R1" :
                            this->isSecondMate() ? 
                                "_R2" :
                                "_R?") :
                qName;
    }
    
    bool isSplicedRead() {
        for(CigarOp op : cigar) {
            if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
                return true;
            }
        }
        
        return false;
    }
    
    uint32_t getNbJunctionsInRead() {
        
        int32_t nbJunctions = 0;
        for(CigarOp op : cigar) {
            if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
                nbJunctions++;
            }
        }
        
        return nbJunctions;
    }
    
    bool isMultiplySplicedRead() {        
        return getNbJunctionsInRead() > 1;
    }
    
    uint32_t calcNbAlignedBases() {
        
        uint32_t count = 0;
        for (CigarOp op : cigar) {
            
            switch ( op.type ) {

                // for 'M', 'I', '=', 'X' - write bases
                case (BAM_CIGAR_MATCH_CHAR) :
                case (BAM_CIGAR_INS_CHAR) :
                case (BAM_CIGAR_EQUAL_CHAR) :
                case (BAM_CIGAR_DIFF_CHAR) :
                case (BAM_CIGAR_DEL_CHAR) :
                case (BAM_CIGAR_PAD_CHAR) :
                case (BAM_CIGAR_REFSKIP_CHAR) :
                    count += op.length;                    
                    break;
            }
        }
        
        return count;
    }
    
    uint16_t calcMinimalMatchInCigarDataSubset(uint32_t start, uint32_t end) {
        
        if (start > position + calcNbAlignedBases() || end < position)
            BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                    "Found an alignment that does not have a presence in the requested region")));
        
        uint16_t mismatches = 0;
        
        for(CigarOp op : cigar) {
           
            if (position > end) {
                break;
            }
            
            if (CigarOp::opConsumesReference(op.type)) {
                position += op.length;
            }
            
            if (position >= start && op.type == BAM_CIGAR_DIFF_CHAR) {
                mismatches++;
            }
        } 
        
        return end - start - mismatches;
    }
    
    uint16_t alignedBasesBetween(int32_t start, int32_t end) {
        
        if (start > position + calcNbAlignedBases() || end < position)
            BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                    "Found an alignment that does not have a presence in the requested region")));
        
        int32_t length = 0;
        uint32_t pos = position;
        for(CigarOp op : cigar) {
           
            if (pos > end) {
                break;
            }
            
            if (CigarOp::opConsumesReference(op.type)) {
                pos += op.length;
                if (pos >= start && op.type != BAM_CIGAR_REFSKIP_CHAR) {
                    length += op.length;
                }
            }            
        } 
        
        return length;
    }
    
};

typedef shared_ptr<BamAlignment> BamAlignmentPtr;

struct RefSeq {
    int32_t index;
    string name;
    uint32_t length;
    
    RefSeq() : RefSeq(-1, "", 0) {};
    
    RefSeq(const int32_t _index, const string& _name, const uint32_t _length) {
        index = _index;
        name = _name;
        length = _length;
    }
};

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
    
    
public:
    
    BamReader(const path& _bamFile, uint16_t threads);
    
    virtual ~BamReader();
    
    vector<RefSeq> getRefs() const { return refs; }
    
    bam_hdr_t* getHeader() const { return header; }
    
    void open();
    
    void close();
    
    bool next();
    
    /**
     * Returns an unmanaged BamAlignment pointer.  By unmanaged this means that
     * the contents of the pointer will be overwritten when "next()" is called.
     * Because we return a shared pointer the BamAlignment wrapper will also be
     * deleted when it goes out of scope (assuming no other references are added).
     * @return Shared pointer to current bam alignment
     */
    BamAlignmentPtr current();
        
    void setRegion(const int32_t seqIndex, const int32_t start, const int32_t end);
    
    bool isCoordSortedBam();
};


class BamWriter {
private:
    path bamFile;
    
    BGZF *fp;    

public:
    BamWriter(const path& _bamFile) {
        bamFile = _bamFile;
        
    }
    
    virtual ~BamWriter() {}
    
    void open(bam_hdr_t* header) {
        // split
        fp = bam_open(bamFile.c_str(), "w");
        if (fp == NULL) {
            BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                    "Could not open output BAM file: ") + bamFile.string()));
        }
        
        if (bam_hdr_write(fp, header) != 0) {
            BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                    "Could not write header into: ") + bamFile.string()));
        }
    }
    
    int write(BamAlignmentPtr ba) {
         return bam_write1(fp, ba->getRaw());       
    }
    
    void close() {
        bam_close(fp);
    }
};
    
// ***** Static stuff
   
class SamtoolsHelper {
        
public:
    
    static path samtoolsExe;
    
    static bool isCoordSortedBam(const path& bamFile);
    
    /**
     * Creates a command that can be used to merge multiple BAM files with samtools
     * @param bamFiles The paths to each BAM file to merge
     * @param mergedBamFile The output file
     * @param threads Number of threads to use during merging
     * @return Command line
     */
    static string createMergeBamCmd(const vector<path>& bamFiles, const path& mergedBamFile, uint16_t threads);
    
    
    /**
     * Creates a samtools command that can be used to sort a bam file.  Assumes 
     * entries will be sorted by position using 1 thread and 1GB or RAM.
     * @param unsortedFile The bam file that needs sorting
     * @param sortedFile The path to the new sorted bam file which will be created     
     * @return The command that can be used to sort the bam file
     */
    static string createSortBamCmd(const path& unsortedFile, const path& sortedFile) {
        return createSortBamCmd(unsortedFile, sortedFile, false, 1, string("1G"));
    }
    
    /**
     * Creates a samtools command that can be used to sort a bam file
     * @param unsortedFile The bam file that needs sorting
     * @param sortedFile The path to the new sorted bam file which will be created
     * @param sortByName If true, bam entries are sorted by name, otherwise by position
     * @param threads Number of threads to use
     * @param memory Amount of memory to request
     * @return The command that can be used to sort the bam file
     */
    static string createSortBamCmd( const path& unsortedFile, 
                                    const path& sortedFile, 
                                    bool sortByName, 
                                    uint16_t threads, 
                                    const string& memory);
    /**
     * Creates a samtools command that can be used to index a sorted bam file
     * @param sortedBam Path to a sorted bam file to index
     * @return The command that can be used to index the sorted bam file
     */
    static string createIndexBamCmd(const path& sortedBam);
    
};
}
