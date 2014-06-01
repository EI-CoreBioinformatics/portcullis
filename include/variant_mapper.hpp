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

#include <bcf.h>

namespace portculis {

typedef boost::error_info<struct VariantMapperError,string> VariantMapperErrorInfo;
struct VariantMapperException: virtual boost::exception, virtual std::exception { };
    
class VariantMapper {
private:
 
    // Path to the original variant file in bcf format
    string bcfFile;
    bool forcePrep;
    bool verbose;
    
    // Handle to genome map.  Created by constructor.
    bcf_idx_t* bcfIndex;
    
protected:
    
    
public:
    
    /**
     * Creates a Variant Mapper object for a genome file in fasta format.  This
     * uses Samtools to create a fasta index for the genome file and then
     * manages the data structure returned after loading the index.
     */
    VariantMapper(string _bcfFile, bool _forcePrep, bool _verbose) : 
        bcfFile(_bcfFile), forcePrep(_forcePrep), verbose(_verbose) {
            bcfIndex = NULL;
    }
    
    virtual ~VariantMapper() {        
        
        if (bcfIndex != NULL) {
            bcf_idx_destroy(bcfIndex);
        }
    }
        
    string getBcfIndexFile() const {
        return bcfFile + ".bci";
    }
    
    
    void buildBcfIndex() {
        int bcfRes = bcf_idx_build(bcfFile.c_str());

        if (bcfRes != 0) {
            BOOST_THROW_EXCEPTION(GenomeMapperException() << GenomeMapperErrorInfo(string(
                    "BCF indexing failed: ") + bcfFile));
        }
    }
    
        
    void loadBcfIndex() {
        
        if (!exists(bcfFile)) {
           BOOST_THROW_EXCEPTION(GenomeMapperException() << GenomeMapperErrorInfo(string(
                    "Genome file does not exist: ") + bcfFile)); 
        }
        
        string bcfIndexFile = getBcfIndexFile();
        if (!exists(bcfIndexFile)) {
           BOOST_THROW_EXCEPTION(GenomeMapperException() << GenomeMapperErrorInfo(string(
                    "BCF index file does not exist: ") + bcfIndexFile)); 
        }
        
        // Load indices
        bcfIndex = bcf_idx_load(bcfFile.c_str());
    }
    
    
    uint64_t bcfQuery(int tid, int beg) {
        return bcf_idx_query(bcfIndex, tid, beg);
    }
        
};
}

