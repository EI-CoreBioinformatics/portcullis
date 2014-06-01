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

#include <string>
#include <vector>
using std::string;
using std::vector;

#include <boost/exception/all.hpp>

#include <api/BamReader.h>;
using namespace BamTools;

namespace portculis {

typedef boost::error_info<struct DepthParserError,string> DepthParserErrorInfo;
struct DepthParserException: virtual boost::exception, virtual std::exception { };
    
struct DepthLine {
    string ref;
    int32_t pos;
    uint32_t depth;
};

class DepthParser {
private:
 
    // Path to the original genome file in fasta format
    string depthFile;
    
    ifstream* ifs;
    RefVector* refs;
    
    DepthLine last;
    bool start;
    
    
protected:
   
    DepthLine parseLine(string& line) {
        if ( !line.empty() ) {
            vector<string> parts; // #2: Search for tokens
            boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

            if (parts.size() != 3) {
                BOOST_THROW_EXCEPTION(DepthParserException() << DepthParserErrorInfo(string(
                    "Malformed depth file: ") + depthFile));
            }

            return DepthLine(parts[0], lexical_cast<int32_t>(parts[1]), lexical_cast<uint32_t>(parts[2]));
        }
    }
    
    size_t findSize(string& refName) {
        
    }
    
public:
    
    DepthParser(string _depthFile, RefVector* _refs) : depthFile(_depthFile), refs(_refs) {
    
        ifs = new ifstream(depthFile.c_str());
        start = true;
    }
    
    virtual ~DepthParser() {
        
        if (ifs != NULL) {
            delete ifs;
        }
    }
    
    bool loadNextBatch(vector<uint32_t>& depths) {
        
        depths.clear();
        
        string line;
        string thisRef;
        string lastRef;
        if (first) {
            std::getline(ifs, line);
            
            DepthLine dl = parseLine(line);
            
            last = dl;
        }
        
        // Create the vector
        depths.resize(findSize(last.ref), 0);
        
        // Use the details from the last run
        depths[last.pos] = last.depth;
        
        
        // Loop through until end of file or we move onto the next ref seq
        while ( std::getline(ifs, line) ) {
            
            DepthLine dl = parseLine(line);
            
            if (last.ref == dl.ref) {
                depths[dl.pos] = dl.depth;
            }
            else {
                last = dl;
                break;
            }
        }
    }
    
};
}