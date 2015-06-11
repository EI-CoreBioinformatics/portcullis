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
#include <bam.h>

#include "samtools_helper.hpp"
using portcullis::BamAlignment;
using portcullis::BamAlignmentPtr;


// ****** BamReader methods *********

portcullis::BamReader::BamReader(const path& _bamFile, const uint16_t _threads) {
    bamFile = _bamFile;
    threads = _threads;
    
    header = nullptr;    
    index = nullptr;
    iter = nullptr;
    
    c = nullptr;    
}

portcullis::BamReader::~BamReader() {
    
    if (header != nullptr) {
        bam_hdr_destroy(header);
    }
    
    if (c != nullptr) {
        bam_destroy1(c);
    }
        
    if (index != nullptr) {
        hts_idx_destroy(index);
    }
    
    if (iter != nullptr) {
        free(iter->off); 
        free(iter->bins.a); 
        free(iter);
    }
}

void portcullis::BamReader::open() {
    
    // split
    fp = bam_open(bamFile.c_str(), "r");
    if (fp == NULL) {
        BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Could not open input BAM files: ") + bamFile.string()));
    }
    
    // Load header
    header = bam_hdr_read(fp);    
    
    // Identify all the reference sequences in the BAM
    for(uint32_t i = 0; i < header->n_targets; i++) {
        refs.push_back(RefSeq(i, string(header->target_name[i]), header->target_len[i]));
    }
    
    // Load the index
    index = bam_index_load(bamFile.c_str());
    
    // Initialise an empty bam alignment
    c = bam_init1();
}

void portcullis::BamReader::close() {
    bam_close(fp);
}

bool portcullis::BamReader::next(BamAlignment& al) {
    bool res = bam_iter_read(fp, iter, c) >= 0;
    al.setRaw(c);    
    return res;    
}

BamAlignmentPtr portcullis::BamReader::current() {
    return make_shared<BamAlignment>(c, false);
}

void portcullis::BamReader::setRegion(const int32_t seqIndex, const int32_t start, const int32_t end) {
    iter = sam_itr_queryi(index, seqIndex, start, end);
}
    

bool portcullis::BamReader::isCoordSortedBam() {    
    string headerText = header->text;   
    return headerText.find("SO:coordinate") != std::string::npos;
}


// ******* Depth parser methods ********

int portcullis::DepthParser::read_bam(void *data, bam1_t *b) {
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
    if (!(b->core.flag&BAM_FUNMAP)) {
        if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
        else if (aux->min_len && bam_cigar2qlen(&b->core, bam1_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
    }

    return ret;
}

int portcullis::DepthParser::read_bam_skip_gapped(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret = 0;
    bool skip = false;
    do {
        skip = false;
        ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
        uint32_t *cigar = bam1_cigar(b);            
        for(int k=0; k < b->core.n_cigar; ++k) {
            int cop = cigar[k] & BAM_CIGAR_MASK; // operation
            if (cop == BAM_CREF_SKIP) {
                skip = true;
                break;
            }
        }

        if (!(b->core.flag&BAM_FUNMAP)) {
            if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
            else if (aux->min_len && bam_cigar2qlen(&b->core, cigar) < aux->min_len) b->core.flag |= BAM_FUNMAP;
        }

    } while(skip);
    return ret;
}
    
    
portcullis::DepthParser::DepthParser(path _bamFile, uint8_t _strandSpecific, bool _allowGappedAlignments) : 
    bamFile(_bamFile), strandSpecific(_strandSpecific), allowGappedAlignments(_allowGappedAlignments) {

    data = (aux_t**)calloc(1, sizeof(aux_t**));
    data[0] = (aux_t*)calloc(1, sizeof(aux_t));
    data[0]->fp = bam_open(bamFile.c_str(), "r");
    data[0]->min_mapQ = 0;                    
    data[0]->min_len  = 0;                 

    header = bam_header_read(data[0]->fp);
    mplp = allowGappedAlignments ? 
        bam_mplp_init(1, read_bam, (void**)data) :
        bam_mplp_init(1, read_bam_skip_gapped, (void**)data);

    res = 0;
    start = true;
}
    
portcullis::DepthParser::~DepthParser() {

    bam_mplp_destroy(mplp);

    bam_header_destroy(header);

    bam_close(data[0]->fp);
    if (data[0]->iter) { 
        bam_iter_destroy(data[0]->iter);
    }
    free(data);
}
    
         

bool portcullis::DepthParser::loadNextBatch(vector<uint32_t>& depths) {

    if (res == 0 && !start) {
        return false;
    }

    depths.clear();

    int pos = 0;
    int tid = -1;
    int n_plp = 0; // n_plp is the number of covering reads from the i-th BAM

    // the core multi-pileup loop
    const bam_pileup1_t** plp = (const bam_pileup1_t**)calloc(1, sizeof(void*)); // plp points to the array of covering reads (internal in mplp)

    if (start) {
        start = false;

        if ((res = bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp)) > 0) {

            int m = 0;
            for (int j = 0; j < n_plp; ++j) {
                const bam_pileup1_t *p = plp[0] + j;
                if (p->is_del || p->is_refskip) ++m;
            }
            last.ref = tid;
            last.pos = pos+1;
            last.depth = n_plp - m; 
        }
        else {
            free(plp);
            return false;
        }
    }

    // Create the vector
    depths.resize(header->target_len[last.ref], 0);

    // Use the details from the last run
    depths[last.pos] = last.depth;

    while ((res = bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp)) > 0) {

        int m = 0;
        for (int j = 0; j < n_plp; ++j) {
            const bam_pileup1_t *p = plp[0] + j;
            if (p->is_del || p->is_refskip) ++m;
        }

        int32_t rpos = pos+1;
        uint32_t cnt = n_plp - m;

        if (last.ref == tid) {
            // Set the depth
            depths[rpos] = cnt;
        }
        else {

            last.ref = tid;
            last.pos = rpos;
            last.depth = cnt;
            break;
        }
    }

    free(plp);

    return true;
}


// ******** Genome mapper ********

/**
 * Creates a Genome Mapper object for a genome file in fasta format.  This
 * uses Samtools to create a fasta index for the genome file and then
 * manages the data structure returned after loading the index.
 */
portcullis::GenomeMapper::GenomeMapper(path _genomeFile) : 
    genomeFile(_genomeFile) {
        fastaIndex = nullptr;
}

portcullis::GenomeMapper::~GenomeMapper() {        

    if (fastaIndex != nullptr) {
        fai_destroy(fastaIndex);
    }
}
    
    
    
/**
 * Constructs the index for this fasta genome file
 */
void portcullis::GenomeMapper::buildFastaIndex() {

    int faiRes = fai_build(genomeFile.c_str());

    if (faiRes != 0) {
        BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Genome indexing failed: ") + genomeFile.string()));
    }
}
    
/**
 * Loads the index for this genome file.  This must be done before using any
 * of the fetch commands.
 */
void portcullis::GenomeMapper::loadFastaIndex() {

    if (!exists(genomeFile)) {
       BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Genome file does not exist: ") + genomeFile.string())); 
    }

    path fastaIndexFile = getFastaIndexFile();
    if (!exists(fastaIndexFile)) {
       BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Genome index file does not exist: ") + fastaIndexFile.string())); 
    }

    fastaIndex = fai_load(genomeFile.c_str());
}
 



// ****** Samtools Helper methods *********

path portcullis::SamtoolsHelper::samtoolsExe = "samtools";

bool portcullis::SamtoolsHelper::isCoordSortedBam(const path& bamFile) {

    BGZF *fp;
    
    // split
    fp = strcmp(bamFile.c_str(), "-")? bgzf_open(bamFile.c_str(), "r") : bgzf_dopen(fileno(stdin), "r");
    if (fp == NULL) {
        BOOST_THROW_EXCEPTION(SamtoolsException() << SamtoolsErrorInfo(string(
                "Could not open input BAM files: ") + bamFile.string()));
    }
    
    bam_hdr_t* header = bam_hdr_read(fp);
    
    string headerText = header->text;
   
    bool found = false;
    if (headerText.find("SO:coordinate") != std::string::npos) {
        found = true;        
    }
    
    bam_hdr_destroy(header);
        
    return found;
}

/**
 * Creates a command that can be used to merge multiple BAM files with samtools
 * @param samtoolsExe The path to samtools
 * @param bamFiles The paths to each BAM file to merge
 * @param mergedBamFile The output file
 * @param threads Number of threads to use during merging
 * @return Command line
 */
string portcullis::SamtoolsHelper::createMergeBamCmd(const vector<path>& bamFiles, 
                                                     const path& mergedBamFile, 
                                                     uint16_t threads) {

    stringstream inputFiles;
    for(path p : bamFiles) {            
        inputFiles << " " << p.c_str();
    }

    return samtoolsExe.string() + " merge -f -@ " + lexical_cast<string>(threads) + 
            " " + mergedBamFile.string() + 
            inputFiles.str();
}
    
    
/**
 * Creates a samtools command that can be used to sort a bam file
 * @param samtoolsExe The path to samtools
 * @param unsortedFile The bam file that needs sorting
 * @param sortedFile The path to the new sorted bam file which will be created
 * @param sortByName If true, bam entries are sorted by name, otherwise by position
 * @param threads Number of threads to use
 * @param memory Amount of memory to request
 * @return The command that can be used to sort the bam file
 */
string portcullis::SamtoolsHelper::createSortBamCmd(const path& unsortedFile, 
                                                    const path& sortedFile, 
                                                    bool sortByName, 
                                                    uint16_t threads, 
                                                    const string& memory) {

    return samtoolsExe.string() + " sort -@ " + lexical_cast<string>(threads) + 
            " -m " + memory + " " + (sortByName ? "-n " : "") + unsortedFile.string() + 
            " " + sortedFile.string();
}
    
/**
 * Creates a samtools command that can be used to index a sorted bam file
 * @param samtoolsExe The path to samtools
 * @param sortedBam Path to a sorted bam file to index
 * @return The command that can be used to index the sorted bam file
 */
string portcullis::SamtoolsHelper::createIndexBamCmd(const path& sortedBam) {

    return samtoolsExe.string() + " index " + sortedBam.string();        
}



    