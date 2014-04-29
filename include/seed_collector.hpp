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

#include <vector>

#include <boost/foreach.hpp>

#include <api/BamReader.h>
#include <api/BamMultiReader.h>

using std::vector;

using namespace BamTools;

namespace portculis {

class SeedCollector {
private:
 
    vector<string> bamFiles;    
    bool collectiveMode;
    
    vector<BamReader> singleReaders;
    BamMultiReader multiReader;
    
    SamHeader header;
    RefVector refs;
    
protected:
    
    
public:
    
    SeedCollector(vector<string> _bamFiles, bool _collectiveMode) : 
            bamFiles(_bamFiles), collectiveMode(_collectiveMode) {
        
        if (collectiveMode) {
            if (!multiReader.Open(bamFiles)) {
                throw "Could not open bam files";
            }
            
            header = reader.GetHeader();
            refs = reader.GetReferenceData();
        }
        else {
            BOOST_FOREACH(string bamFile, bamFiles) {
                BamReader reader;                
                if (!reader.Open(bamFile)) {
                    throw "Could not open bam files";
                }
                singleReaders.push_back(reader);
            }
        }                
        
    }
            
    void sort() {
        
    }
};
}

