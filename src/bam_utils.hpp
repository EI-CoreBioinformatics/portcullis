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

#include <sstream>
#include <string>
#include <vector>
using namespace::std;

#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
using boost::filesystem::exists;
using boost::lexical_cast;

#include <api/algorithms/Sort.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;

namespace portcullis {
namespace bamtools {

#ifdef SAMTOOLS_PATH
    const string SAMTOOLS_EXE = string(SAMTOOLS_PATH) + "/samtools";
#else
    const string SAMTOOLS_EXE("samtools");
#endif    
    
typedef boost::error_info<struct BamError,string> BamErrorInfo;
struct BamException: virtual boost::exception, virtual std::exception { };

class BamUtils {
public:
    
    static void mergeBams(vector<string>& bamFiles, string bamOut) {
        
        BamMultiReader reader;
        if (!reader.Open(bamFiles)) {
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo("Could not open input BAM files"));
        }

        const SamHeader header = reader.GetHeader();
        const RefVector refs = reader.GetReferenceData();
        
        BamWriter writer;
        if (!writer.Open(bamOut, header, refs)) {
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string("Could not open output BAM file for merged results: ") + bamOut));
        }
        
        BamAlignment al;
        while (reader.GetNextAlignment(al)) {
            writer.SaveAlignment(al);
        }
        
        reader.Close();
        writer.Close();
    }
    
    static bool isSortedBam(string bamFile) {
        
        BamReader reader;
        
        if (!reader.Open(bamFile)) {
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string("Could not open input BAM files: ") + bamFile));
        }

        const SamHeader header = reader.GetHeader();
        
        reader.Close();
        
        return header.HasSortOrder();
    }
    
    /**
     * Creates a samtools command that can be used to sort a bam file.  Assumes 
     * entries will be sorted by position using 1 thread and 1GB or RAM.
     * @param unsortedFile The bam file that needs sorting
     * @param sortedFile The path to the new sorted bam file which will be created     
     * @return The command that can be used to sort the bam file
     */
    static string createSortBamCmd(string unsortedFile, string sortedFile) {
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
    static string createSortBamCmd(string unsortedFile, string sortedFile, bool sortByName, uint16_t threads, string memory) {
        
        return SAMTOOLS_EXE + " sort -@ " + lexical_cast<string>(threads) + 
                " -m " + memory + " " + (sortByName ? "-n " : "") + unsortedFile + 
                " " + sortedFile;
    }
    
    /**
     * Creates a samtools command that can be used to index a sorted bam file
     * @param sortedBam Path to a sorted bam file to index
     * @return The command that can be used to index the sorted bam file
     */
    static string createIndexBamCmd(string sortedBam) {
        
        return SAMTOOLS_EXE + " index " + sortedBam;        
    }
    
    static bool opAlignsToReference(char type) {
        
        switch (type) {
            // increase end position on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                return true;
            default:
                return false;
        }
    }
    
    static bool opFollowsReference(char type) {
        
        switch (type) {
            // increase end position on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                return true;
            default:
                return false;
        }
    }
    
    static string deriveName(const BamAlignment& al) {
        
        return al.IsPaired() ? 
                al.Name + (al.IsFirstMate() ? 
                            "_R1" :
                            al.IsSecondMate() ? 
                                "_R2" :
                                "_R?") :
                al.Name;
    }
    
    static uint16_t calcMinimalMatchInCigarDataSubset(BamAlignment& ba, int32_t start, int32_t end) {
        
        if (start > ba.Position + ba.AlignedBases.size() || end < ba.Position)
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                    "Found an alignment that does not have a presence in the requested region")));
        
        int32_t pos = ba.Position;
        uint16_t mismatches = 0;
        int32_t length = 0;
        bool inRegion = false;
        
        for(CigarOp op : ba.CigarData) {
           
            if (pos > end) {
                break;
            }
            
            if (BamUtils::opFollowsReference(op.Type)) {
                pos += op.Length;
            }
            
            if (pos >= start && op.Type == 'X') {
                mismatches++;
            }
        } 
        
        return end - start - mismatches;
    }
    
    static uint16_t alignedBasesBetween(BamAlignment& ba, int32_t start, int32_t end) {
        
        if (start > ba.Position + ba.AlignedBases.size() || end < ba.Position)
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                    "Found an alignment that does not have a presence in the requested region")));
        
        int32_t pos = ba.Position;
        int32_t length = 0;
        
        for(CigarOp op : ba.CigarData) {
           
            if (pos > end) {
                break;
            }
            
            if (BamUtils::opFollowsReference(op.Type)) {
                pos += op.Length;
                if (pos >= start && op.Type != Constants::BAM_CIGAR_REFSKIP_CHAR) {
                    length += op.Length;
                }
            }            
        } 
        
        return length;
    }
};
}    
}
