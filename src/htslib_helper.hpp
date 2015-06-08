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

typedef boost::error_info<struct HtslibError,string> HtslibErrorInfo;
struct HtslibException: virtual boost::exception, virtual std::exception { };

class HtslibHelper {
   
public:

    
    
    
        
    
    
    
    
    static bool isSplicedRead(BamAlignment& ba) {
        for(CigarOp op : ba.CigarData) {
            if (op.Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
                return true;
            }
        }
        
        return false;
    }
    
    static int32_t getNbJunctionsInRead(BamAlignment& ba) {
        
        int32_t nbJunctions = 0;
        for(CigarOp op : ba.CigarData) {
            if (op.Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
                nbJunctions++;
            }
        }
        
        return nbJunctions;
    }
    
    static bool isMultiplySplicedRead(BamAlignment& ba) {
        
        return getNbJunctionsInRead(ba) > 1;
    }
    
    static uint16_t calcMinimalMatchInCigarDataSubset(BamAlignment& ba, int32_t start, int32_t end) {
        
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
    
    static uint16_t alignedBasesBetween(BamAlignment& ba, int32_t start, int32_t end) {
        
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
};
}    
}
