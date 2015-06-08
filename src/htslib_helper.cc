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

#include <sstream>
#include <string>
#include <vector>
using std::vector;
using std::stringstream;
using std::string;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::lexical_cast;

#include <bam.h>

#include "htslib_helper.hpp"

    

    

bool portcullis::HtslibHelper::opAlignsToReference(char type) {

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

bool portcullis::HtslibHelper::opFollowsReference(char type) {

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

string portcullis::HtslibHelper::deriveName(const BamAlignment& al) {

    return al.IsPaired() ? 
            al.Name + (al.IsFirstMate() ? 
                        "_R1" :
                        al.IsSecondMate() ? 
                            "_R2" :
                            "_R?") :
            al.Name;
}

bool portcullis::HtslibHelper::isSplicedRead(BamAlignment& ba) {
    for(CigarOp op : ba.CigarData) {
        if (op.Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
            return true;
        }
    }

    return false;
}

int32_t portcullis::HtslibHelper::getNbJunctionsInRead(BamAlignment& ba) {

    int32_t nbJunctions = 0;
    for(CigarOp op : ba.CigarData) {
        if (op.Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
            nbJunctions++;
        }
    }

    return nbJunctions;
}

bool portcullis::HtslibHelper::isMultiplySplicedRead(BamAlignment& ba) {

    return getNbJunctionsInRead(ba) > 1;
}

uint16_t portcullis::HtslibHelper::calcMinimalMatchInCigarDataSubset(BamAlignment& ba, int32_t start, int32_t end) {

    if (start > ba.Position + ba.AlignedBases.size() || end < ba.Position)
        BOOST_THROW_EXCEPTION(BamUtilsException() << BamUtilsErrorInfo(string(
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

uint16_t portcullis::HtslibHelper::alignedBasesBetween(BamAlignment& ba, int32_t start, int32_t end) {

    if (start > ba.Position + ba.AlignedBases.size() || end < ba.Position)
        BOOST_THROW_EXCEPTION(BamUtilsException() << BamUtilsErrorInfo(string(
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
