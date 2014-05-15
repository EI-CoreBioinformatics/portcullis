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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "portculis.hpp"
#include "genome_mapper.hpp"
#include "seed_collector.hpp"
#include "bam_utils.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;


using portculis::GenomeMapper;
using portculis::SeedCollector;
using portculis::sortBam;
using portculis::mergeBams;
using portculis::indexBam;

namespace po = boost::program_options;


// Default values for arguments
const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;
const uint32_t DEFAULT_CHUNK_SIZE_PER_THREAD = 10000;
const uint32_t DEFAULT_GAP_SIZE = 100;


std::pair<string, string> bamPrep(vector<string> bamFiles, string outputPrefix, bool forcePrep) {
    
    string indexedBamFile;
    string sortedBamFile;
    if (bamFiles.size() > 1) {
        string mergedBam = outputPrefix + ".merged.bam";
        string sortedBam = outputPrefix + ".merged.sorted.bam";
        string indexedBam = outputPrefix + ".merged.sorted.bam.bti";

        bool mergedBamExists = boost::filesystem::exists(mergedBam);
        
        if (mergedBamExists) {
            cout << "Pre-merged BAM detected: " << mergedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing merge anyway due to user request." << endl;
            }
        }
        
        if (!mergedBamExists || forcePrep) {
            cout << "Found " << bamFiles.size() << " BAM files." << endl 
                 << "Merging BAMS...";
            mergeBams(bamFiles, mergedBam);
            cout << "done." << endl;
        }

        bool mergedAndSortedBamExists = boost::filesystem::exists(sortedBam);
        
        if (mergedAndSortedBamExists) {
            cout << "Pre-merged and sorted BAM detected: " << sortedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing sort anyway due to user request." << endl;
            }
        }
        
        if (!mergedAndSortedBamExists || forcePrep) {
            cout << "Sorting " << mergedBam << " ... ";
            sortBam(mergedBam, sortedBam, false);
            cout << "done." << endl;        
        }

        bool indexedBamExists = boost::filesystem::exists(indexedBam);
        
        if (indexedBamExists) {
            cout << "Pre-merged and sorted BAM index detected: " << indexedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing index creation anyway due to user request." << endl;
            }
        }
        
        if (!indexedBamExists || forcePrep) {
            cout << "Indexing " << sortedBam << " ... ";
            indexBam(sortedBam, indexedBam);
            cout << "done." << endl;
        }

        indexedBamFile = indexedBam;
        sortedBamFile = sortedBam;
    }
    else {  // 1 bam File

        string bamFile = bamFiles[0];
        string bamFileToIndex = bamFile;

        // Sort the BAM file if necessary
        bool bamIsSorted = portculis::isSortedBam(bamFile);

        if (bamIsSorted) {
            cout << "Pre-sorted BAM detected: " << bamFile << endl;
            
            if (forcePrep) {
                cout << "Forcing sort anyway due to user request." << endl;
            }
        }
        
        if (!bamIsSorted || forcePrep) {
            string sortedBam = bamFile + ".sorted.bam";

            cout << "Sorting " << bamFile << " ... ";
            sortBam(bamFile, sortedBam, false);
            cout << "done" << endl;

            bamFileToIndex = sortedBam;
        }

        string indexedBam = bamFile + ".sorted.bam.bti";

        bool bamIsIndexed = boost::filesystem::exists(indexedBam);

        if (bamIsIndexed) {
            cout << "Pre-indexed BAM detected: " << indexedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing index creation anyway due to user request." << endl;
            }
        }

        // Index the BAM file if necessary.
        if (!bamIsIndexed || forcePrep) {
            cout << "Indexing " << bamFileToIndex << " ... ";
            indexBam(bamFileToIndex, indexedBam);
            cout << "done." << endl;
        }
        indexedBamFile = indexedBam;
        sortedBamFile = bamFileToIndex;
    }
    
    return std::pair<string, string>(sortedBamFile, indexedBamFile);
}


/**
 * Start point for portculis.
 */
int main(int argc, char *argv[]) {
    try {
        // Portculis args
        vector<string> bamFiles;
        string genomeFile;
        string outputPrefix;
        //bool collectiveMode;
        bool forcePrep;
        uint16_t threads;
        bool disableThreadedIO;
        uint32_t gap_size;
        bool verbose;
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options("Portculis Help.\nUsage: portculis [options] <genome-file> (<bam-file>)+ \nAllowed options");
        generic_options.add_options()
                ("output_prefix,o", po::value<string>(&outputPrefix)->default_value(DEFAULT_OUTPUT_PREFIX), "Path prefix for files generated by this program.")
                //("collective,c", po::value<bool>(&collectiveMode)->default_value(false), "Whether to treat all provided bam files as one large merged bam file or to handle each separately.") 
                ("force_prep,f", po::bool_switch(&forcePrep)->default_value(false), "Whether to force preparation (sorting and indexing) of the input bam files.  By default, portculis only sorts bam files if the header does not contain a sorted (coordinate) flag.  Also, it only indexes if it cannot find a bam file with a bamtools bti extension suffix.")
                ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_THREADS), "The number of threads to use.")
                ("disable_threaded_io", po::value<bool>(&disableThreadedIO)->default_value(false), "Whether to acquire data in parallel to processing.")
                ("gap,g", po::value<uint32_t>(&gap_size)->default_value(DEFAULT_GAP_SIZE), "The minimum gap size between adjacent alignments that determines whether or not we can chunk the data at this point.")
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print extra information")
                ("version", po::bool_switch(&version)->default_value(false), "Print version string")
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

        // Output version information then exit if requested
        if (version) {
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "Portculis"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.0"
#endif
            cout << PACKAGE_NAME << " V" << PACKAGE_VERSION << endl;
            return 0;
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

        // Print configuration information if in verbose mode
        if (verbose) {
            
            cerr << "Portculis settings: " << endl
                    << " Provided " << bamFiles.size() << " bam files" << endl
                    << " Output prefix: " << outputPrefix << endl
                    << " Worker threads: " << threads << endl
                    << " Use threaded IO: " << !disableThreadedIO << endl
                    << " Forcing BAM prep: " << forcePrep << endl
                    << " Minimum gap size: " << gap_size << endl;
        }

        // OK, we're good to do some real work now!

        // Create a map of the genome (wrapper for faidx from samtools)
        GenomeMapper genomeMapper(genomeFile);
        
        // Prep the BAM input to produce a usable sorted bam plus bamtools bti index
        std::pair<string, string> bamPrepResult = bamPrep(bamFiles, outputPrefix, forcePrep);
        
        // First job is to collect seeds (i.e. alignments containing an 'N' in their cigar)
        string sortedBam = bamPrepResult.first;
        string indexedBam = bamPrepResult.second;
        string seedFile = outputPrefix + ".seeds.bam";
        
        cout << "Collecting seeds from " << sortedBam << " ... ";
        uint64_t seedCount = SeedCollector(sortedBam, seedFile, verbose).collect();
        cout << "done. Found " << seedCount << " seeds." << endl;
        
        if (seedCount > 0) {
        
            //Portculis instance(bam_file, genome_file, output_prefix, threads, !disable_threaded_io, chunk_size_per_thread, gap_size, verbose);
            portculis::Portculis instance;
            instance.process();
        }
        else {
            
            // No point in carrying on if no seeds were found
            cout << "No seeds found.  Portculis finishing." << endl;
        }

    } catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 4;
    } catch (const char* msg) {
        cerr << "Error: " << msg << endl;
        return 5;
    } catch (...) {
        cerr << "Error: Exception of unknown type!" << endl;
        return 6;
    }

    return 0;
}

