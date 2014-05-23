//  ********************************************************************
//  This file is part of Portculis.
//
//  Portculis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portculis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <sstream>
#include <string>
#include <vector>
using namespace::std;

#include <boost/exception/all.hpp>
#include <boost/foreach.hpp>

#include <api/algorithms/Sort.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;

#include <bamtools_sort.h>
       

namespace portculis {
namespace bamtools {

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
    
    static void sortBam(string unsortedFile, string sortedFile, bool sortByName) {
        
        SortSettings settings;
        settings.inputBamFilename = unsortedFile;
        settings.outputBamFilename = sortedFile;
        settings.isSortingByName = sortByName;
        
        SortTool sorter(&settings);
        
        if (!sorter.run()) {
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                    "Failed to successfully sort: ") + unsortedFile));
        }  
    }
    
    static void indexBam(string sortedBam) {
        
        BamReader reader;
        
        if (!reader.Open(sortedBam)) {
            BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                    "Could not open BAM file to index: ") + sortedBam));
        }
        
        reader.CreateIndex(BamIndex::BAMTOOLS);        
        reader.Close();
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
};
}    
}
