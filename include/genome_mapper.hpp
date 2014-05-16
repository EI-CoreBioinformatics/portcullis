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
#include <iostream>

#include <boost/filesystem.hpp>

#include <faidx.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::string;

namespace portculis {

class GenomeMapper {
private:
 
    // Path to the original genome file in fasta format
    string genomeFile;
    bool forcePrep;
    
    // Handle to genome map.  Created by constructor.
    faidx_t* index;
    
protected:
    
    
public:
    
    /**
     * Creates a Genome Mapper object for a genome file in fasta format.  This
     * uses Samtools to create a fasta index for the genome file and then
     * manages the data structure returned after loading the index.
     */
    GenomeMapper(string _genomeFile, bool _forcePrep) : 
        genomeFile(_genomeFile), forcePrep(_forcePrep) {
        
        string indexFile = genomeFile + string(".fai");
                
        bool indexExists = boost::filesystem::exists(indexFile);
        
        if (indexExists) {
            cout << "Indexed genome detected: " << indexFile << endl;
            
            if (forcePrep) {
                cout << "Forcing index build anyway due to user request." << endl;
            }
        }
        
        if (!indexExists || forcePrep) {
            int res = fai_build(genomeFile.c_str());

            if (res != 0)
                throw "Could not build the fasta index for the genome file.";
        }
        
        index = fai_load(genomeFile.c_str());
        
        cout << "Genome index acquired" << endl;
    }
    
    virtual ~GenomeMapper() {        
        fai_destroy(index);
    }
    
    /**
     * @abstract    Fetch the sequence in a region.
     * @param  reg  Region in the format "chr2:20,000-30,000"
     * @param  len  Length of the region
     * @return      Pointer to the sequence; null on failure
     * @discussion The returned sequence is allocated by malloc family
     * and should be destroyed by end users by calling free() on it.
     */
    char* fetch(const char* reg, int* len) {
        return fai_fetch(index, reg, len);        
    }
    
    /**
     * @abstract    Fetch the sequence in a region.
     * @param  name Region name
     * @param  start    Start location on region (zero-based)
     * @param  end  End position (zero-based)
     * @param  len  Length of the region
     * @return      Pointer to the sequence; null on failure
     * @discussion The returned sequence is allocated by malloc family
     * and should be destroyed by end users by calling free() on it.
     */
    char* fetch(char* name, int start, int end, int* len) {
        return faidx_fetch_seq(index, name, start, end, len);        
    }
    
    int getNbSeqs() {
        return faidx_fetch_nseq(index); 
    }
    
    string getFaidxPath() {
        return genomeFile + ".fai";
    }
    
};
}

