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

#include <glob.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using boost::algorithm::trim;
namespace bfs = boost::filesystem;
using bfs::path;
namespace po = boost::program_options;

#include "bam_utils.hpp"
#include "genome_mapper.hpp"
#include "portcullis_fs.hpp"
using portcullis::bamtools::BamUtils;
using portcullis::PortcullisFS;


namespace portcullis {
    
typedef boost::error_info<struct PrepareError,string> PrepareErrorInfo;
struct PrepareException: virtual boost::exception, virtual std::exception { };


const string DEFAULT_PREP_OUTPUT_DIR = "portcullis_prep_data";
const uint16_t DEFAULT_PREP_THREADS = 1;

const string PORTCULLIS = "portcullis";

const string FASTA_EXTENSION = ".fa";
const string FASTA_INDEX_EXTENSION = ".fai";
const string BAM_EXTENSION = ".bam";
const string BAM_INDEX_EXTENSION = ".bai";
const string BCF_EXTENSION = ".bcf";
const string BCF_INDEX_EXTENSION = ".bci";
const string BAM_DEPTH_EXTENSION = ".bdp";
   
enum class StrandSpecific : std::uint8_t {
    UNSTRANDED = 0,
    FIRSTSTRAND = 1,
    SECONDSTRAND = 2
};

inline const string SSToString(StrandSpecific ss) {
    switch (ss) {
        case StrandSpecific::UNSTRANDED:   return "UNSTRANDED";
        case StrandSpecific::FIRSTSTRAND:  return "FIRSTSTRAND";
        case StrandSpecific::SECONDSTRAND: return "SECONDSTRAND";
        default:      return "[Unknown StrandSpecific type]";
    }
}

inline const StrandSpecific SSFromString(string& ss) {
    
    if (boost::iequals(ss, "UNSTRANDED")) {
        return StrandSpecific::UNSTRANDED;
    }
    else if (boost::iequals(ss, "FIRSTSTRAND")) {
        return StrandSpecific::FIRSTSTRAND;
    }
    else if (boost::iequals(ss, "SECONDSTRAND")) {
        return StrandSpecific::SECONDSTRAND;
    }
        
    BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Can't recognise StrandSpecific string: ") + ss));
    
    return StrandSpecific::UNSTRANDED;
}

class PreparedFiles {
    
private:
    path prepDir;
    
public:
    
    PreparedFiles(const path& _prepDir) : prepDir(_prepDir) {
        
        if (!bfs::exists(prepDir)) {
            bfs::create_directory(prepDir);
        }
    }
        
    path getPrepDir() const {
        return prepDir;
    }

    path getUnsortedBamFilePath() const {
        return path(prepDir.string() + "/" + PORTCULLIS + ".unsorted.alignments" + BAM_EXTENSION);
    }
    
    path getSortedBamFilePath() const {
        return path(prepDir.string() + "/" + PORTCULLIS + ".sorted.alignments" + BAM_EXTENSION);
    }
    
    path getBamIndexFilePath() const {
        return path(getSortedBamFilePath().string() + BAM_INDEX_EXTENSION);
    }
    
    path getBcfFilePath() const {
        return path(getSortedBamFilePath().string() + BCF_EXTENSION);
    }    
    
    path getBcfIndexFilePath() const {
        return path(getBcfFilePath().string() + BCF_INDEX_EXTENSION);
    }    
    
    path getGenomeFilePath() const {
        return path(prepDir.string() + "/" + PORTCULLIS + ".genome" + FASTA_EXTENSION);
    }
    
    path getGenomeIndexFilePath() const {
        return path(getGenomeFilePath().string() + FASTA_INDEX_EXTENSION);
    }
    
    path getSettingsFilePath() const {
        return path(prepDir.string() + "/" + PORTCULLIS + ".settings");
    }
    
    bool valid();
    
    void clean();
    
    StrandSpecific loadSettings();
};


class Prepare {

private:
    
    PreparedFiles* output;
    StrandSpecific strandSpecific;
    bool force;
    bool useLinks;
    uint16_t threads;
    bool verbose;
    PortcullisFS fs;

    
public:
    
    Prepare(const path& _outputPrefix) : Prepare(_outputPrefix, StrandSpecific::UNSTRANDED, false, false, 1, false) {
    }
    
    Prepare(const path& _outputPrefix, StrandSpecific _strandSpecific, bool _force, bool _useLinks, uint16_t _threads, bool _verbose);
    
    virtual ~Prepare() {
        delete output;
    }
    
    
protected:
    
    bool copy(const path& from, const path& to, const string& msg);
    
    bool genomeIndex();
    
    
    /**
     * Merge together a set of BAM files, use the output prefix to construct a
     * file name for the merged file
     * @param bamFiles The set of bam files to merge
     * @return 
     */
    bool bamMerge(vector<path> bamFiles);

    /**
     * Sorts the unsorted bam file if required or forced
     * @param inputBam
     * @return 
     */
    bool bamSort();
    
    bool bamIndex();

public:
    
    static vector<path> globFiles(vector<path> input);
    
    void prepare(vector<path> bamFiles, const path& originalGenomeFile);
    
    PortcullisFS getFs() const {
        return fs;
    }

    void setFs(PortcullisFS fs) {
        this->fs = fs;
    }

    

    bool outputDetails();

  
    static string helpMessage() {
        return string("\nPortcullis Prepare Mode Help.\n\n") +
                      "Usage: portcullis prep [options] <genome-file> (<bam-file>)+ \n\n" +
                      "Allowed options";
    }
    
    
    
    static int main(int argc, char *argv[], PortcullisFS& fs);
    
};
}