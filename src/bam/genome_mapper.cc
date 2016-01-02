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

#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::make_shared;
using std::shared_ptr;
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
#include <htslib/sam.h>

#include "genome_mapper.hpp"

// ******** Genome mapper ********

/**
 * Creates a Genome Mapper object for a genome file in fasta format.  This
 * uses Samtools to create a fasta index for the genome file and then
 * manages the data structure returned after loading the index.
 */
portcullis::bam::GenomeMapper::GenomeMapper(path _genomeFile) : 
    genomeFile(_genomeFile) {
        fastaIndex = nullptr;
}

portcullis::bam::GenomeMapper::~GenomeMapper() {        

    if (fastaIndex != nullptr) {
        fai_destroy(fastaIndex);
    }
}
    
    
    
/**
 * Constructs the index for this fasta genome file
 */
void portcullis::bam::GenomeMapper::buildFastaIndex() {

    int faiRes = fai_build(genomeFile.c_str());

    if (faiRes != 0) {
        BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Genome indexing failed: ") + genomeFile.string()));
    }
}
    
/**
 * Loads the index for this genome file.  This must be done before using any
 * of the fetch commands.
 */
void portcullis::bam::GenomeMapper::loadFastaIndex() {

    if (!exists(genomeFile)) {
       BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Genome file does not exist: ") + genomeFile.string())); 
    }

    path fastaIndexFile = getFastaIndexFile();
    if (!exists(fastaIndexFile)) {
       BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
                "Genome index file does not exist: ") + fastaIndexFile.string())); 
    }

    fastaIndex = fai_load(genomeFile.c_str());
}
 

/**
* @abstract    Fetch the sequence in a region.
* @param  reg  Region in the format "chr2:20,000-30,000"
* @param  len  Length of the region
* @return      The sequence as a string; empty string if no seq found
*/
string portcullis::bam::GenomeMapper::fetchBases(const char* reg, int* len) const {
   char* cseq = fai_fetch(fastaIndex, reg, len);
   string strseq = cseq == NULL ? string("") : string(cseq);
   if (cseq != NULL)
       free(cseq);
   return strseq;        
}

/**
 * @abstract    Fetch the sequence in a region.
 * @param  name Region name
 * @param  start    Start location on region (zero-based, inclusive)
 * @param  end  End position (zero-based, inclusive)
 * @param  len  Length of the region
 * @return      The sequence as a string; empty string if no seq found
 */
string portcullis::bam::GenomeMapper::fetchBases(const char* name, int start, int end, int* len) const {
    char* cseq = faidx_fetch_seq(fastaIndex, name, start, end, len);
    string strseq = cseq == NULL ? string("") : string(cseq);
    if (cseq != NULL)
        free(cseq);
    return strseq;        
}