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
using std::vector;
using std::stringstream;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::lexical_cast;

namespace portcullis {

typedef boost::error_info<struct SamtoolsError,string> SamtoolsErrorInfo;
struct SamtoolsException: virtual boost::exception, virtual std::exception { };

class SamtoolsHelper {
        
public:
    
    // CIGAR constants
    static const uint8_t BAM_CIGAR_MATCH    = 0;
    static const uint8_t BAM_CIGAR_INS      = 1;
    static const uint8_t BAM_CIGAR_DEL      = 2;
    static const uint8_t BAM_CIGAR_REFSKIP  = 3;
    static const uint8_t BAM_CIGAR_SOFTCLIP = 4;
    static const uint8_t BAM_CIGAR_HARDCLIP = 5;
    static const uint8_t BAM_CIGAR_PAD      = 6;
    static const uint8_t BAM_CIGAR_SEQMATCH = 7;
    static const uint8_t BAM_CIGAR_MISMATCH = 8;

    static const char BAM_CIGAR_MATCH_CHAR    = 'M';
    static const char BAM_CIGAR_INS_CHAR      = 'I';
    static const char BAM_CIGAR_DEL_CHAR      = 'D';
    static const char BAM_CIGAR_REFSKIP_CHAR  = 'N';
    static const char BAM_CIGAR_SOFTCLIP_CHAR = 'S';
    static const char BAM_CIGAR_HARDCLIP_CHAR = 'H';
    static const char BAM_CIGAR_PAD_CHAR      = 'P';
    static const char BAM_CIGAR_SEQMATCH_CHAR = '=';
    static const char BAM_CIGAR_MISMATCH_CHAR = 'X';    
    
    
private:
    
    path bamFile;
    
    
// ***** Static stuff
    
public:
    
    
    /**
     * Creates a command that can be used to merge multiple BAM files with samtools
     * @param samtoolsExe The path to samtools
     * @param bamFiles The paths to each BAM file to merge
     * @param mergedBamFile The output file
     * @param threads Number of threads to use during merging
     * @return Command line
     */
    static string createMergeBamCmd(const path& samtoolsExe, const vector<path>& bamFiles, const path& mergedBamFile, uint16_t threads);
    
    
    /**
     * Creates a samtools command that can be used to sort a bam file.  Assumes 
     * entries will be sorted by position using 1 thread and 1GB or RAM.
     * @param samtoolsExe The path to samtools
     * @param unsortedFile The bam file that needs sorting
     * @param sortedFile The path to the new sorted bam file which will be created     
     * @return The command that can be used to sort the bam file
     */
    static string createSortBamCmd(const path& samtoolsExe, const path& unsortedFile, const path& sortedFile) {
        return createSortBamCmd(samtoolsExe, unsortedFile, sortedFile, false, 1, string("1G"));
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
    static string createSortBamCmd( const path& samtoolsExe, 
                                    const path& unsortedFile, 
                                    const path& sortedFile, 
                                    bool sortByName, 
                                    uint16_t threads, 
                                    string memory);
    /**
     * Creates a samtools command that can be used to index a sorted bam file
     * @param samtoolsExe The path to samtools
     * @param sortedBam Path to a sorted bam file to index
     * @return The command that can be used to index the sorted bam file
     */
    static string createIndexBamCmd(const path& samtoolsExe, const path& sortedBam);
    
    
    
    static bool isCoordSortedBam(const path& bamFile);
    
    static inline bool opAlignsToReference(char type) {
        
        switch (type) {
            // increase end position on CIGAR chars [DMXN=]
            case BAM_CIGAR_MATCH_CHAR    :
            case BAM_CIGAR_MISMATCH_CHAR :
            case BAM_CIGAR_SEQMATCH_CHAR :
                return true;
            default:
                return false;
        }
    }
    
    static inline bool opFollowsReference(char type) {
        
        switch (type) {
            // increase end position on CIGAR chars [DMXN=]
            case BAM_CIGAR_DEL_CHAR      :
            case BAM_CIGAR_MATCH_CHAR    :
            case BAM_CIGAR_MISMATCH_CHAR :
            case BAM_CIGAR_REFSKIP_CHAR  :
            case BAM_CIGAR_SEQMATCH_CHAR :
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
    
};
}
