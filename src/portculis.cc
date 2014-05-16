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
#include "bam_prepare.hpp"
#include "seed_collector.hpp"


using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;

using portculis::Portculis;
using portculis::GenomeMapper;
using portculis::SeedCollector;


namespace po = boost::program_options;


// Default values for arguments
const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;
const uint32_t DEFAULT_CHUNK_SIZE_PER_THREAD = 10000;
const uint32_t DEFAULT_GAP_SIZE = 100;


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
        cout << endl 
             << "Indexing genome" << endl
             << "---------------" << endl;
        GenomeMapper genomeMapper(genomeFile, forcePrep);
        
        // Prep the BAM input to produce a usable sorted bam plus bamtools bti index
        cout << endl
             << "Preparing BAM(s)" << endl
             << "----------------" << endl;
        std::pair<string, string> bamPrepResult = portculis::bamPrep(bamFiles, outputPrefix, forcePrep);
        
        // First job is to collect seeds (i.e. alignments containing an 'N' in their cigar)
        string sortedBam = bamPrepResult.first;
        string indexedBam = bamPrepResult.second;
        string seedFile = outputPrefix + string(".seeds.bam");
        string unsplicedFile = outputPrefix + string(".unspliced.bam");
        
        cout << endl 
             << "Seed collecting" << endl
             << "---------------" << endl;
        portculis::CollectorResults results;
        SeedCollector(sortedBam, seedFile, unsplicedFile, verbose).collect(results);
        results.report(cout);
        
        if (results.seedCount > 0) {
        
            cout << endl 
                 << "Processing seeds" << endl
                 << "----------------" << endl;
            Portculis instance(seedFile, sortedBam, genomeFile, outputPrefix, threads, verbose);
            instance.process();
            cout << endl << "Portculis finished" << endl;
        }
        else {            
            // No point in carrying on if no seeds were found
            cout << "WARNING: No seeds found." << endl
                 << "Deleting unspliced alignments file to save disk space (unspliced file will be identical to sorted input)." << endl
                 << "Portculis finishing." << endl;
            
            if (!boost::filesystem::remove(unsplicedFile)) {
                throw "Error deleting unspliced file";
            }
        }
        
        cout.flush();

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

