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

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::string;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using boost::filesystem::absolute;
using boost::filesystem::copy_file;
using boost::filesystem::remove;
using boost::filesystem::exists;
using boost::filesystem::create_symlink;
using boost::filesystem::create_directory;
using boost::filesystem::symbolic_link_exists;

#include "bam_utils.hpp"
#include "genome_mapper.hpp"
using portculis::bamtools::BamUtils;


namespace portculis {
    
typedef boost::error_info<struct PrepareError,string> PrepareErrorInfo;
struct PrepareException: virtual boost::exception, virtual std::exception { };


const string DEFAULT_PREP_OUTPUT_DIR = "portculis_prep_data";
const uint16_t DEFAULT_PREP_THREADS = 1;

const string PORTCULIS = "portculis";

const string FASTA_EXTENSION = ".fa";
const string FASTA_INDEX_EXTENSION = ".fai";
const string BAM_EXTENSION = ".bam";
const string BAM_INDEX_EXTENSION = ".bti";
const string BCF_EXTENSION = ".bcf";
const string BCF_INDEX_EXTENSION = ".bci";
const string BAM_DEPTH_EXTENSION = ".bdp";
    

class PreparedFiles {
    
private:
    string prepDir;
    
public:
    
    PreparedFiles(string _prepDir) : prepDir(_prepDir) {
        
        if (!exists(prepDir)) {
            create_directory(prepDir);
        }
    }
        
    string getPrepDir() const {
        return prepDir;
    }

    string getUnsortedBamFilePath() const {
        return prepDir + "/" + PORTCULIS + ".unsorted.alignments" + BAM_EXTENSION;
    }
    
    string getSortedBamFilePath() const {
        return prepDir + "/" + PORTCULIS + ".sorted.alignments" + BAM_EXTENSION;
    }
    
    string getBamIndexFilePath() const {
        return getSortedBamFilePath() + BAM_INDEX_EXTENSION;
    }
    
    string getBamDepthFilePath() const {
        return getSortedBamFilePath() + BAM_DEPTH_EXTENSION;
    }

    string getBcfFilePath() const {
        return getSortedBamFilePath() + BCF_EXTENSION;
    }    
    
    string getBcfIndexFilePath() const {
        return getBcfFilePath() + BCF_INDEX_EXTENSION;
    }    
    
    string getGenomeFilePath() const {
        return prepDir + "/" + PORTCULIS + ".genome" + FASTA_EXTENSION;
    }
    
    string getGenomeIndexFilePath() const {
        return getGenomeFilePath() + FASTA_INDEX_EXTENSION;
    }
    
    string getSettingsFilePath() const {
        return prepDir + "/" + PORTCULIS + ".settings";
    }
    
    bool valid() {
        
        if (!exists(getSortedBamFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find sorted BAM files at: ") + getSortedBamFilePath()));
        }

        if (!exists(getBamIndexFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find BAM index at: ") + getBamIndexFilePath()));
        }
        
        if (!exists(getBamDepthFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find BAM depth file at: ") + getBamDepthFilePath()));
        }
        
        if (!exists(getGenomeFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find genome file at: ") + getGenomeFilePath()));
        }
        
        if (!exists(getGenomeIndexFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find genome index at: ") + getGenomeIndexFilePath()));
        }
        
        if (!exists(getSettingsFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find settings file at: ") + getSettingsFilePath()));
        } 
        
        return true;
    }
    
    void clean() {
        
        remove(getUnsortedBamFilePath());
        remove(getSortedBamFilePath());
        remove(getBamIndexFilePath());
        remove(getBamDepthFilePath());
        remove(getGenomeFilePath());
        remove(getGenomeIndexFilePath());
        remove(getBcfFilePath());
        remove(getBcfIndexFilePath());
        remove(getSettingsFilePath());
    }
};


class Prepare {

private:
    
    PreparedFiles* output;
    bool strandSpecific;
    bool force;
    bool useLinks;
    uint16_t threads;
    bool verbose;

    void init(string _outputDir, bool _strandSpecific, bool _force, bool _useLinks, uint16_t _threads, bool _verbose) {
        
        output = new PreparedFiles(_outputDir);
        strandSpecific = _strandSpecific;
        force = _force;
        useLinks = _useLinks;
        threads = _threads;
        verbose = _verbose;
        
        if (verbose) {
            cout << "Configured portculis prep to use the following settings: " << endl
                 << " - Output directory: " << output->getPrepDir() << endl
                 << " - Strand specific library: " << boolalpha << strandSpecific << endl
                 << " - Force prep (cleans output directory): " << boolalpha << force << endl
                 << " - Use symbolic links instead of copy where possible: " << boolalpha << useLinks << endl
                 << " - Threads (for sorting BAM): " << threads << endl << endl;
        }
        
        if (force) {
            if (verbose) {
                cout << "Cleaning output dir " << _outputDir << " ... ";
                cout.flush();
            }
            output->clean();
            if (verbose) cout << "done." << endl << endl;
        }
    }
    
public:
    
    Prepare(string _outputPrefix) {
        init(_outputPrefix, false, false, false, 1, false);
    }
    
    Prepare(string _outputPrefix, bool _strandSpecific, bool _force, bool _useLinks, uint16_t _threads, bool _verbose) {
        init(_outputPrefix, _strandSpecific, _force, _useLinks, _threads, _verbose);
    }
    
    virtual ~Prepare() {
        delete output;
    }
    
    
protected:
    
    bool copy(string from, string to, string msg) {
        
        auto_cpu_timer timer(1, string(" - Copy ") + msg + " - Wall time taken: %ws\n\n");
        
        bool fileExists = exists(to);
        
        if (fileExists) {
            if (verbose) cout << "Prepped " << msg << " file detected: " << to << endl;            
        }
        else {
            
            if (useLinks) {
                create_symlink(absolute(from), to);            
                if (verbose) cout << "Created symlink from " << from << " to " << to << endl;
            }
            else {
                if (verbose) {
                    cout << "Copying from " << from << " to " << to << " ... ";
                    cout.flush();
                }
                copy_file(from, to);
                if (verbose) cout << "done." << endl;
            }
        }
        
        return exists(to) || symbolic_link_exists(to);
    }
    
    bool genomeIndex() {
        
        auto_cpu_timer timer(1, " - Genome Index - Wall time taken: %ws\n\n");
        
        const string genomeFile = output->getGenomeFilePath();
        const string bcfFile = output->getBcfFilePath();
        const string indexFile = output->getGenomeIndexFilePath();
        
        bool indexExists = exists(indexFile);
        
        if (indexExists) {
            if (verbose) cout << "Pre-indexed genome detected: " << indexFile << endl;            
        }
        else {
            
            if (verbose) {
                cout << "Indexing genome " << genomeFile << " ... ";
                cout.flush();
            }
            
            // Create the index
            GenomeMapper(genomeFile, bcfFile, force, verbose).buildFastaIndex();
            
            if (verbose) {
                cout << "done." << endl
                     << "Genome index file created at: " << output->getGenomeIndexFilePath() << endl;                
            }
        }
        
        return exists(indexFile);
    }
    
    
    /**
     * Merge together a set of BAM files, use the output prefix to construct a
     * file name for the merged file
     * @param bamFiles The set of bam files to merge
     * @return 
     */
    bool bamMerge(vector<string> bamFiles) {

        auto_cpu_timer timer(1, " - BAM Merge - Wall time taken: %ws\n\n");
        
        string mergedBam = output->getUnsortedBamFilePath();        

        bool mergedBamExists = exists(mergedBam);

        if (mergedBamExists) {
            if (verbose) cout << "Pre-merged BAM detected: " << mergedBam << endl;
        }
        else {
            
            if (verbose) {
                cout << "Found " << bamFiles.size() << " BAM files." << endl 
                     << "Merging BAMS...";
                cout.flush();
            }
            
            BamUtils::mergeBams(bamFiles, mergedBam);
            
            if (verbose) cout << "done." << endl;
        }
        
        // Return true if the merged BAM exists now, which is should do
        return exists(mergedBam);        
    }

    /**
     * Sorts the unsorted bam file if required or forced
     * @param inputBam
     * @return 
     */
    bool bamSort() {
        
        auto_cpu_timer timer(1, " - BAM Sort - Wall time taken: %ws\n\n");
        
        const string unsortedBam = output->getUnsortedBamFilePath();
        const string sortedBam = output->getSortedBamFilePath();
        
        bool sortedBamExists = exists(sortedBam);
        
        if (sortedBamExists) {            
            if (verbose) cout << "Pre-sorted BAM detected: " << sortedBam << endl;
        }
        else {
            
            if (BamUtils::isSortedBam(unsortedBam) && !force) {
                
                if (verbose) cout << "Provided BAM appears to be sorted already, just creating symlink instead." << endl;
                create_symlink(unsortedBam, sortedBam);            
                if (verbose) cout << "Created symlink from " << unsortedBam << " to " << sortedBam << endl;
            }
            else {
                
                if (verbose) {
                    cout << "Sorting " << unsortedBam << " ... " << endl;
                }

                // Sort the BAM file by coordinate
                BamUtils::sortBam(unsortedBam, sortedBam, false, threads, "1G");
            }
            
        }
        
        // Return true if the sorted BAM exists now, which is should do
        return exists(sortedBam);        
    }
    
    bool bamIndex() {
        
        auto_cpu_timer timer(1, " - BAM Index - Wall time taken: %ws\n\n");
        
        const string sortedBam = output->getSortedBamFilePath();
        const string indexedFile = output->getBamIndexFilePath();
        
        bool indexedBamExists = exists(indexedFile);
        
        if (indexedBamExists) {
            if (verbose) cout << "Pre-indexed BAM detected: " << indexedFile << endl;
        }
        else {
            if (verbose) {
                cout << "Indexing " << sortedBam << " ... ";
                cout.flush();
            }
            
            // Create BAM index
            BamUtils::indexBam(sortedBam);
            
            if (verbose) {
                cout << "done." << endl
                     << "BAM index created at: " << output->getBamIndexFilePath() << endl;
            }
        }
        
        return exists(indexedFile);
    }
    
    bool bamDepth() {
        
        auto_cpu_timer timer(1, " - BAM Depth - Wall time taken: %ws\n\n");
        
        const string sortedBam = output->getSortedBamFilePath();
        const string depthFile = output->getBamDepthFilePath();
        
        if (exists(depthFile)) {
            if (verbose) cout << "Depth file detected: " << depthFile << endl;            
        }
        else {
            
            if (verbose) {
                cout << "Calculating per base depth from " << sortedBam << " ... " << endl;
            }
            
            // Create BAM pileup
            BamUtils::depth(sortedBam, depthFile);
            
        }
        
        return exists(depthFile);
    }
    
    bool bamPileup() {
        
        auto_cpu_timer timer(1, " - BAM Pileup - Wall time taken: %ws\n\n");
        
        const string sortedBam = output->getSortedBamFilePath();
        const string genomeFile = output->getGenomeFilePath();
        string pileupFile = output->getBcfFilePath();
        
        bool pileupFileExists = exists(pileupFile);
        
        if (pileupFileExists) {
            if (verbose) cout << "Pre-piled-up file detected: " << pileupFile << endl;            
        }
        else {
            
            if (verbose) {
                cout << "Piling up alignments from " << sortedBam << " ... " << endl;
            }
            
            // Create BAM pileup
            BamUtils::pileupBam(sortedBam, genomeFile, pileupFile);
            
        }
        
        return exists(pileupFile);
    }

    bool bcfIndex() {
        
        auto_cpu_timer timer(1, " - BCF Index - Wall time taken: %ws\n\n");
        
        const string bcfPileups = output->getBcfFilePath();
        const string bcfIndex = output->getBcfIndexFilePath();
        const string genomeFile = output->getGenomeFilePath();
        
        bool bcfIndexExists = exists(bcfIndex);
        
        if (bcfIndexExists) {
            if (verbose) cout << "BCF pileup index file detected: " << bcfIndex << endl;            
        }
        else {
            
            if (verbose) {
                cout << "Indexing BCF pileup file " << bcfPileups << " ... " << endl;
            }
            
            // Create BCF index
            GenomeMapper(genomeFile, bcfPileups, force, verbose).buildBcfIndex();
            
            if (verbose) {
                cout << "done." << endl
                     << "BCF index created at: " << bcfIndex << endl;
            }
        }
        
        return exists(bcfIndex);
    }

public:
    
    
    void prepare(vector<string> bamFiles, string originalGenomeFile) {

        const bool doMerge = bamFiles.size() > 1;        
        
        if (bamFiles.empty()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "No BAM files to process")));
        }

        string mergedBamFile = doMerge ? output->getUnsortedBamFilePath() : bamFiles[0];
        
        // Merge the bams to output unsorted bam file if required, otherwise just
        // copy / symlink the file provided
        if (doMerge) {
            if (!bamMerge(bamFiles)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not merge BAM files")));
            }
        }
        else {
            
            // Copy / Symlink the file to the output dir
            if (!copy(bamFiles[0], output->getUnsortedBamFilePath(), "BAM")) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                            "Could not copy/symlink BAM file to: ") + output->getUnsortedBamFilePath()));
            }            
        }

        // Sort the file
        if (!bamSort()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not sort: ") + output->getUnsortedBamFilePath()));
        }

        // Index the sorted file
        if (!bamIndex()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Failed to index: ") + output->getSortedBamFilePath()));
        }
        
        // Test if depth exists
        if (!bamDepth()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not create depth file")));
        }
        
        // Copy / Symlink the file to the output dir
        if (!copy(originalGenomeFile, output->getGenomeFilePath(), "genome")) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not copy/symlink genome file to: ") + output->getGenomeFilePath()));
        }
        
        // Test if index exists
        if (!genomeIndex()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not create genome map")));
        }
        
        /*// Create pileups
        if (!bamPileup()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not pileup: ") + output->getSortedBamPath()));
        }
        
        // Create pileup index
        if (!bcfIndex()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not create index for: ") + output->getBcfFilePath()));
        }*/
    }

    bool outputDetails() {
        
        const string settingsFile = output->getSettingsFilePath();
        
        std::ofstream outfile(settingsFile.c_str());
        if(!outfile.is_open()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not open settings file for writing: ") + settingsFile));
        }

        outfile << "SS=" << (strandSpecific ? "true" : "false") << endl;
        outfile.close();
    }

  
    static string helpMessage() {
        return string("\nPortculis Prepare Mode Help.\n\n") +
                      "Usage: portculis prep [options] <genome-file> (<bam-file>)+ \n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portculis args
        vector<string> bamFiles;
        string genomeFile;
        string outputDir;
        bool force;
        bool strandSpecific;
        bool useLinks;
        uint16_t threads;
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output,o", po::value<string>(&outputDir)->default_value(DEFAULT_PREP_OUTPUT_DIR), 
                    (string("Output directory for prepared files. Default: ") + DEFAULT_PREP_OUTPUT_DIR).c_str())
                ("force,f", po::bool_switch(&force)->default_value(false), 
                    "Whether or not to clean the output directory before processing, thereby forcing full preparation of the genome and bam files.  By default portculis will only do what it thinks it needs to.")
                ("strand_specific,ss", po::bool_switch(&strandSpecific)->default_value(false), 
                    "Whether BAM alignments were generated using a strand specific RNAseq library.")
                ("use_links,l", po::bool_switch(&useLinks)->default_value(false), 
                    "Whether to use symbolic links from input data to prepared data where possible.  Saves time and disk space but is less robust.")
                ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_PREP_THREADS),
                    (string("The number of threads to used to sort the BAM file (if required).  Default: ") + lexical_cast<string>(DEFAULT_PREP_THREADS)).c_str())
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                    "Print extra information")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("bam-files,i", po::value< vector<string> >(&bamFiles), "Path to the BAM files to process.")
                ("genome-file,g", po::value<string>(&genomeFile), "Path to the genome file to process.")
                ;

        // Positional option for the input bam file
        po::positional_options_description p;
        p.add("genome-file", 1);
        p.add("bam-files", 100);

        // Combine non-positional options
        po::options_description cmdline_options;
        cmdline_options.add(generic_options).add(hidden_options);

        // Parse command line
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);

        // Output help information the exit if requested
        if (help || argc <= 1) {
            cout << generic_options << endl;
            return 1;
        }
        
        // Acquire path to bam file
        if (vm.count("bam-files")) {
            bamFiles = vm["bam-files"].as<vector<string> >();
        }

        // Acquire path to genome file
        if (vm.count("genome-file")) {
            genomeFile = vm["genome-file"].as<string>();
        }
        
        // Test if provided genome exists
        if (!exists(genomeFile)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not find genome file at: ") + genomeFile));
        }

        auto_cpu_timer timer(1, "\nPortculis prep completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portculis in prepare mode" << endl
             << "---------------------------------" << endl << endl;
        
        // Create the prepare class
        Prepare prep(outputDir, strandSpecific, force, useLinks, threads, verbose);
        
        // Prep the input to produce a usable indexed and sorted bam plus, indexed
        // genome and queryable coverage information
        prep.prepare(bamFiles, genomeFile);
        
        // Output any remaining details
        prep.outputDetails();
        
    }
};
}