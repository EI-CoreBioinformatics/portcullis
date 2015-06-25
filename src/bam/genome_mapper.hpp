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

#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::shared_ptr;
using std::make_shared;
using std::string;
using std::vector;
using std::stringstream;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::lexical_cast;

#include <htslib/faidx.h>
#include <bam.h>

#include "bam_master.hpp"

namespace portcullis { namespace bam {

class GenomeMapper {
private:
 
    // Path to the original genome file in fasta format
    path genomeFile;
    
    // Handle to genome map.  Created by constructor.
    faidx_t* fastaIndex;
    
protected:
    
    
public:
    
    /**
     * Creates a Genome Mapper object for a genome file in fasta format.  This
     * uses Samtools to create a fasta index for the genome file and then
     * manages the data structure returned after loading the index.
     */
    GenomeMapper(path _genomeFile);
    
    virtual ~GenomeMapper();
    
    
    path getFastaIndexFile() const {
        return path(genomeFile.parent_path()) /= path(genomeFile.leaf().string() + ".fai");
    }
    
    
    /**
     * Constructs the index for this fasta genome file
     */
    void buildFastaIndex();
    
    /**
     * Loads the index for this genome file.  This must be done before using any
     * of the fetch commands.
     */
    void loadFastaIndex();
    
    
    /**
     * @abstract    Fetch the sequence in a region.
     * @param  reg  Region in the format "chr2:20,000-30,000"
     * @param  len  Length of the region
     * @return      The sequence as a string; empty string if no seq found
     */
    string fetchBases(const char* reg, int* len) const;
    
    /**
     * @abstract    Fetch the sequence in a region.
     * @param  name Region name
     * @param  start    Start location on region (zero-based, inclusive)
     * @param  end  End position (zero-based, exclusive)
     * @param  len  Length of the region
     * @return      The sequence as a string; empty string if no seq found
     */
    string fetchBases(const char* name, int start, int end, int* len) const;
    
    /**
     * Get the number of sequences / contigs / scaffolds in the genome
     * @return 
     */
    int getNbSeqs() {
        return faidx_nseq(fastaIndex); 
    }
    
};

}}
