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
#include <memory>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::shared_ptr;
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

#include <portcullis/bam/bam_master.hpp>
using portcullis::bam::Strandedness;

#include <portcullis/portcullis_fs.hpp>
using portcullis::PortcullisFS;


namespace portcullis {
    
typedef boost::error_info<struct PrepareError,string> PrepareErrorInfo;
struct PrepareException: virtual boost::exception, virtual std::exception { };


const string DEFAULT_PREP_OUTPUT_DIR = "portcullis_prep";
const uint16_t DEFAULT_PREP_THREADS = 1;

const string PORTCULLIS = "portcullis";

const string FASTA_EXTENSION = ".fa";
const string FASTA_INDEX_EXTENSION = ".fai";
const string BAM_EXTENSION = ".bam";
const string BAI_EXTENSION = ".bai";
const string CSI_EXTENSION = ".csi";
const string BCF_EXTENSION = ".bcf";
const string BCF_INDEX_EXTENSION = ".bci";
const string BAM_DEPTH_EXTENSION = ".bdp";
   

class PreparedFiles {
    
private:
    path prepDir;
    
public:
    
    PreparedFiles() {}
    
    PreparedFiles(const path& _prepDir) : prepDir(_prepDir) {    
        
        if (!bfs::exists(prepDir)) {
            if (!bfs::create_directories(prepDir)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not create output directory at: ") + prepDir.string()));
            }
        }
        else if (!bfs::is_directory(prepDir)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "File exists with name of suggested output directory: ") + prepDir.string()));            
        }
    }
    
    void setPrepDir(const path& _prepDir) {
        
        prepDir = _prepDir;
        
        if (!bfs::exists(prepDir)) {
            if (!bfs::create_directories(prepDir)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not create output directory at: ") + prepDir.string()));
            }
        }
        else if (!bfs::is_directory(prepDir)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "File exists with name of suggested output directory: ") + prepDir.string()));            
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
    
    path getBamIndexFilePath(bool useCsi) const {
        return path(getSortedBamFilePath().string() + (useCsi ? CSI_EXTENSION : BAI_EXTENSION));
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
    
    bool valid(bool useCsi) const;
    
    void clean();    
};


class Prepare {

private:
    
    shared_ptr<PreparedFiles> output;
    bool force;
    bool useLinks;
    uint16_t threads;
    bool useCsi;
    bool verbose;
    
    
public:
    
    Prepare(const path& _outputPrefix);
    
    bool isForce() const {
        return force;
    }

    void setForce(bool force) {
        this->force = force;
    }

    shared_ptr<PreparedFiles> getOutput() const {
        return output;
    }

    void setOutput(shared_ptr<PreparedFiles> output) {
        this->output = output;
    }

    uint16_t getThreads() const {
        return threads;
    }

    void setThreads(uint16_t threads) {
        this->threads = threads;
    }

    bool isUseCsi() const {
        return useCsi;
    }

    void setUseCsi(bool useCsi) {
        this->useCsi = useCsi;
    }

    bool isUseLinks() const {
        return useLinks;
    }

    void setUseLinks(bool useLinks) {
        this->useLinks = useLinks;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }

    
    virtual ~Prepare() {       
    }
    
    
protected:
    
    bool copy(const path& from, const path& to, const string& msg, const bool requireFileExists);
    
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
     * @param input Path to input BAM to sort
     * @param output Sorted BAM file
     * @return 
     */
    bool bamSort(const path& input, const path& output);
    
    bool bamIndex(const bool copied);
    
    /**
     * Checks whether the specified indexing method can support the genome sequence lengths
     * @param genomeFile
     * @param useCsi
     * @return 
     */
    bool checkIndexMode(const path& genomeIndexFile, const bool useCsi);

public:
    
    static vector<path> globFiles(vector<path> input);
    
    void prepare(vector<path> bamFiles, const path& originalGenomeFile);
    
    bool outputDetails();

    static string title() {
        return string("Portcullis Prepare Mode Help.");
    }
  
    static string description() {
        return string("Prepares a genome and bam file(s) ready for junction analysis.  This involves\n") +
                      "ensuring the bam file is sorted and indexed and the genome file is indexed.";
    }
    
    static string usage() {
        return string("portcullis prep [options] <genome-file> (<bam-file>)+");
    }
    
    
    static int main(int argc, char *argv[]);
    
};
}