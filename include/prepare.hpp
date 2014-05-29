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
#include <vector>
using std::string;
using std::vector;

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
using boost::timer::auto_cpu_timer;
using boost::filesystem::copy_file;
using boost::filesystem::exists;
using boost::filesystem::create_symlink;

#include "bam_utils.hpp"
using portculis::bamtools::BamUtils;


namespace portculis {
    
typedef boost::error_info<struct PrepareError,string> PrepareErrorInfo;
struct PrepareException: virtual boost::exception, virtual std::exception { };

const string FASTA_EXTENSION = ".fa";
const string FASTA_INDEX_EXTENSION = ".fai";
const string BAM_EXTENSION = ".bam";
const string BAM_INDEX_EXTENSION = ".bti";
    

class PreparedFiles {
    
public:
    static const string PORTCULIS = "portculis";
    
private:
    string prepDir;
    
public:
    
    PreparedFiles(string _prepDir) :
        prepDir(_prepDir) {
    }
        
    string getPrepDir() const {
        return prepDir;
    }

        
    string getBamFile() const {
        return PORTCULIS + ".alignments" + BAM_EXTENSION;
    }
    
    string getBamIndexFile() const {
        return getBamFile() + BAM_INDEX_EXTENSION;
    }
    
    string getGenomeFile() const {
        return PORTCULIS + ".genome" + FASTA_EXTENSION;
    }
    
    string getGenomeIndexFile() const {
        return getGenomeFile() + FASTA_INDEX_EXTENSION;
    }

        
    
};


class Prepare {

public:    
    static const string MERGED_EXTENSION = ".merged";
    static const string SORTED_EXTENSION = ".sorted";
    static const string BAM_EXTENSION = ".bam";
    static const string INDEX_EXTENSION = ".bti";
    static const string PILEUP_EXTENSION = ".pileup";

private:
    
    PreparedFiles* output;
    bool force;
    bool useLinks;
    bool verbose;

    void init(string _outputDir, bool _force, bool _useLinks, bool _verbose) {
        output = new PreparedFiles(_outputDir);
        force = _force;
        useLinks = _useLinks;
        verbose = _verbose;
    }
    
public:
    
    BamPrepare(string _outputPrefix) {
        init(_outputPrefix, false, false, false);
    }
    
    BamPrepare(string _outputPrefix, bool _force, bool _useLinks, bool _verbose) {
        init(_outputPrefix, _force, _useLinks, _verbose);
    }
    
    virtual ~BamPrepare() {
        delete output;
    }
    
    string getMergedBamPath() {
        return output->getPrepDir() + "temp_merged_alignments" + BAM_EXTENSION;
    }
    
    
protected:
    
    /**
     * Merge together a set of BAM files, use the output prefix to construct a
     * file name for the merged file
     * @param bamFiles The set of bam files to merge
     * @return 
     */
    bool merge(vector<string> bamFiles) {

        string mergedBam = getMergedBamPath();        

        bool mergedBamExists = exists(mergedBam);

        if (mergedBamExists && verbose) {
            cout << "Pre-merged BAM detected: " << mergedBam << endl;

            if (force) {
                cout << "Forcing merge anyway due to user request." << endl;
            }
        }

        if (!mergedBamExists || force) {
            
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
     * Sorts the provided bam file if required or forced
     * @param inputBam
     * @return 
     */
    bool sort(string inputBam) {
        
        string sortedBam = inputBam + SORTED_EXTENSION + BAM_EXTENSION;
        
        bool mergedAndSortedBamExists = exists(sortedBam);
        
        if (mergedAndSortedBamExists && verbose) {
            cout << "Pre-merged and sorted BAM detected: " << sortedBam << endl;
            
            if (force) {
                cout << "Forcing sort anyway due to user request." << endl;
            }
        }
        
        if (!mergedAndSortedBamExists || force) {
            if (verbose) {
                cout << "Sorting " << mergedBam << " ... ";
                cout.flush();
            }
            
            // Sort the BAM file by coordinate
            BamUtils::sortBam(inputBam, sortedBam, false);
            
            if (verbose) cout << "done." << endl;        
        }
        
        // Return true if the sorted BAM exists now, which is should do
        return exists(sortedBam);        
    }
    
    bool index(string sortedBam) {
        
        string indexedBam = sortedBam + INDEX_EXTENSION;
        
        bool indexedBamExists = exists(indexedBam);
        
        if (indexedBamExists && verbose) {
            cout << "Pre-merged and sorted BAM index detected: " << indexedBam << endl;
            
            if (force) {
                cout << "Forcing index creation anyway due to user request." << endl;
            }
        }
        
        if (!indexedBamExists || force) {
            
            if (verbose) {
                cout << "Indexing " << sortedBam << " ... ";
                cout.flush();
            }
            
            // Index the BAM
            BamUtils::indexBam(sortedBam);
            
            if (verbose) cout << "done." << endl;
        }
        
        return exists(indexedBam);
    }
    
    bool index(string sortedBam) {
        
        string indexedBam = sortedBam + INDEX_EXTENSION;
        
        bool indexedBamExists = exists(indexedBam);
        
        if (indexedBamExists && verbose) {
            cout << "BAM index detected: " << indexedBam << endl;
            
            if (force) {
                cout << "Forcing indexing anyway due to user request." << endl;
            }
        }
        
        if (!indexedBamExists || force) {
            
            if (verbose) {
                cout << "Indexing " << sortedBam << " ... ";
                cout.flush();
            }
            
            // Index the BAM
            BamUtils::indexBam(sortedBam);
            
            if (verbose) cout << "done." << endl;
        }
        
        return exists(indexedBam);
    }


public:
    
    bool indexGenome(string genomeFile) {
        
        if (useLinks) {
            create_symb
        }
        GenomeMapper genomeMapper(genomeFile);
        
    }
    
    bool prepare(vector<string> bamFiles) {

        const bool doMerge = bamFiles.size() > 1;
        
        auto_cpu_timer timer(1, "Wall time taken: %ws\n");        

        
        if (bamFiles.empty()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "No BAM files to process")));
        }

        string mergedBamFile = doMerge ? getMergedBamPath() : bamFiles[0];
        
        // .. but if there is more than one file then actually do the merging
        if (doMerge) {
            if (!merge(bamFiles)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not merge BAM files")));
            }
        }

        // Sort the file
        if (!sort(mergedBamFile)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not sort: ") + mergedBam));
        }

        // Index the sorted file
        if (!index(output->getBamFile())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Failed to index: ") + mergedBam));
        }
        
        // Create pileups
        if (!pileup(output->getBamFile())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not sort: ") + mergedBam));
        }

    }    

  
    static int main(int argc, char *argv[]) {
        
        // Portculis args
        vector<string> bamFiles;
        string genomeFile;
        string outputDir;
        bool strandSpecific;
        bool useLinks;
        bool verbose;
        bool help;

        // Declare the supported options.
        po::options_description generic_options("Portculis Help.  Prepare Mode.\nUsage: portculis prepare [options] <genome-file> (<bam-file>)+ \nAllowed options");
        generic_options.add_options()
                ("output,o", po::value<string>(&outputDir)->default_value(DEFAULT_OUTPUT_PREFIX), "Output directory for prepared files.")
                ("strand_specific,ss", po::bool_switch(&strandSpecific)->default_value(false), "Whether BAM alignments were generated using a strand specific RNAseq library.")
                ("use_links,l", po::bool_switch(&useLinks)->default_value(false), "Whether to use symbolic links from input data to prepared data where possible.  Saves time and disk space but is less robust.")
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print extra information")
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
        if (help) {
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

        // Test to see if genome file exists
        if (genomeFile == "") {
            cerr << endl << "ERROR: You must specify a genome file!" << endl;
            cerr << endl << generic_options << endl;
            return 2;
        } else if (!boost::filesystem::exists(genomeFile)) {
            cerr << "ERROR: Specified genome file " << genomeFile << " does not exist!" << endl;
            return 2;
        }

        // Test to see if all bam files exists
        if (bamFiles.empty()) {
            cerr << endl << "ERROR: You must specify at least one bam file to process" << endl;
            cout << endl << generic_options << endl;
            return 3;
        } else {

            BOOST_FOREACH(string bamFile, bamFiles) {
                if (!boost::filesystem::exists(bamFile)) {
                    cerr << "ERROR: Specified BAM file " << bamFile << " does not exist!" << endl;
                    return 3;
                }
            }
        }

        // OK, we're good to do some real work now!
        auto_cpu_timer timer(1, "\nTotal runtime: %ws\n");
        
        Prepare prep(outputDir, true, useLinks, verbose);
        
        // Create a map of the genome (wrapper for faidx from samtools)
        cout << endl 
             << "Indexing genome" << endl
             << "---------------" << endl;
        prep.indexGenome(genomeFile);
        
        // Prep the BAM input to produce a usable sorted bam plus bamtools bti index
        cout << endl
             << "Preparing BAM(s)" << endl
             << "----------------" << endl;
        prep.prepare(bamFiles);
    }
};
}