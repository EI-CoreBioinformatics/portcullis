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
using std::cout;
using std::endl;
using std::ifstream;
using std::string;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::exists;
using boost::timer::auto_cpu_timer;

#include <faidx.h>

namespace portculis {

typedef boost::error_info<struct GenomeMapperError,string> GenomeMapperErrorInfo;
struct GenomeMapperException: virtual boost::exception, virtual std::exception { };
    
class GenomeMapper {
private:
 
    // Path to the original genome file in fasta format
    string genomeFile;
    bool forcePrep;
    bool verbose;
    
    // Handle to genome map.  Created by constructor.
    faidx_t* index;
    
protected:
    
    
public:
    
    /**
     * Creates a Genome Mapper object for a genome file in fasta format.  This
     * uses Samtools to create a fasta index for the genome file and then
     * manages the data structure returned after loading the index.
     */
    GenomeMapper(string _genomeFile, bool _forcePrep, bool _verbose) : 
        genomeFile(_genomeFile), forcePrep(_forcePrep), verbose(_verbose) {
    }
    
    virtual ~GenomeMapper() {        
        fai_destroy(index);
    }
    
    
    string getIndexFile() const {
        return genomeFile + ".fai";
    }
    
    /**
     * Constructs the index for this fasta genome file
     */
    void buildIndex() {
        
        int res = fai_build(genomeFile.c_str());

        if (res != 0) {
            BOOST_THROW_EXCEPTION(GenomeMapperException() << GenomeMapperErrorInfo(string(
                    "Genome indexing failed: ") + genomeFile));
        }
    }
    
    /**
     * Loads the index for this genome file.  This must be done before using any
     * of the fetch commands.
     */
    void loadIndex() {
        
        if (!exists(genomeFile)) {
           BOOST_THROW_EXCEPTION(GenomeMapperException() << GenomeMapperErrorInfo(string(
                    "Genome file does not exist: ") + genomeFile)); 
        }
        
        string indexFile = getFaidxPath();
        if (!exists(indexFile)) {
           BOOST_THROW_EXCEPTION(GenomeMapperException() << GenomeMapperErrorInfo(string(
                    "Genome index file does not exist: ") + indexFile)); 
        }
        
        index = fai_load(genomeFile.c_str());
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

