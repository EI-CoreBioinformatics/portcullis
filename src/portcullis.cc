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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
using boost::to_upper_copy;
using boost::timer::auto_cpu_timer;
namespace po = boost::program_options;

#include "bam/bam_master.hpp"
using portcullis::bam::BamHelper;

#include "portcullis_fs.hpp"
#include "junction_builder.hpp"
#include "prepare.hpp"
#include "junction_filter.hpp"
#include "bam_filter.hpp"
#include "cluster.hpp"
#include "train.hpp"
#include "test.hpp"
using portcullis::PortcullisFS;
using portcullis::JunctionBuilder;
using portcullis::Prepare;
using portcullis::JunctionFilter;
using portcullis::BamFilter;
using portcullis::Cluster;
using portcullis::Train;
using portcullis::Test;

typedef boost::error_info<struct PortcullisError,string> PortcullisErrorInfo;
struct PortcullisException: virtual boost::exception, virtual std::exception { };

// Default values for arguments
const uint16_t DEFAULT_THREADS = 4;
const uint32_t DEFAULT_CHUNK_SIZE_PER_THREAD = 10000;
const uint32_t DEFAULT_GAP_SIZE = 100;

enum Mode {
    PREP,
    JUNC,
    FILTER,
    BAM_FILT,
    FULL,
    CLUSTER,
    TRAIN,
    TEST
};

Mode parseMode(string mode) {
    
    string upperMode = boost::to_upper_copy(mode);
    
    if (upperMode == string("PREP")) {
        return PREP;                
    }
    else if (upperMode == string("JUNC")) {
        return JUNC;
    }
    else if (upperMode == string("FILTER")) {
        return FILTER;
    }
    else if (upperMode == string("BAMFILT")) {
        return BAM_FILT;
    }
    else if (upperMode == string("FULL")) {
        return FULL;
    }
    else if (upperMode == string("CLUSTER")) {
        return CLUSTER;
    }
    else if (upperMode == string("TRAIN")) {
        return TRAIN;
    }
    else if (upperMode == string("TEST")) {
        return TEST;
    }
    else {
        BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
                    "Could not recognise mode string: ") + mode));
    }
}

string helpHeader() {
    return string("\nPortcullis Help.\n\n") +
                  "Portcullis is a tool to identify genuine splice junctions using aligned RNAseq reads\n\n" +
                  "Usage: portcullis [options] <mode> <mode_args>\n\n" +
                  "Available modes:\n" +
                  " - prep    - Prepares a genome and bam file(s) ready for junction analysis\n" +
                  " - junc    - Perform junction analysis on prepared data\n" +
                  " - filter  - Discard unlikely junctions and produce BAM containing alignments to genuine junctions\n" +
                  " - bamfilt - Filters a BAM to remove any reads associated with invalid junctions\n" + 
                  " - full    - Runs prep, junc, filter and bamfilt as a complete pipeline\n" +
                  " - cluster - Clusters potential junctions to help distinguish real junctions for false\n\n" +
                  "\nAvailable options";
}

static string fullHelp() {
    return string("\nPortcullis Full Pipeline Mode Help.\n\n") +
                  "Usage: portcullis full [options] <genome-file> (<bam-file>)+ \n\n" +
                  "Allowed options";
}


int mainFull(int argc, char *argv[]) {
    
    // Portcullis args
    std::vector<path> bamFiles;
    path genomeFile;
    path outputDir;
    string strandSpecific;
    uint16_t threads;
    bool useLinks;
    bool force;
    bool useCsi;
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(fullHelp(), 100);
    generic_options.add_options()
            ("output,o", po::value<path>(&outputDir)->default_value("portcullis_out"), 
                "Output directory for prepared files. Default: portcullis_out")
            ("force,f", po::bool_switch(&force)->default_value(false), 
                "Whether or not to clean the output directory before processing, thereby forcing full preparation of the genome and bam files.  By default portcullis will only do what it thinks it needs to.")
            ("strand_specific,ss", po::value<string>(&strandSpecific)->default_value(strandednessToString(Strandedness::UNKNOWN)), 
                "Whether BAM alignments were generated using a strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads)  Default: \"unknown\"")
            ("use_links,l", po::bool_switch(&useLinks)->default_value(false), 
                "Whether to use symbolic links from input data to prepared data where possible.  Saves time and disk space but is less robust.")
            ("use_csi,c", po::bool_switch(&useCsi)->default_value(false), 
                "Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it supports very long target sequences (probably not an issue unless you are working on huge genomes).  BAI has the advantage that it is more widely supported (useful for viewing in genome browsers).")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use.  Default: 1")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("bam-files,i", po::value< std::vector<path> >(&bamFiles), "Path to the BAM files to process.")
            ("genome-file,g", po::value<path>(&genomeFile), "Path to the genome file to process.")
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
        bamFiles = vm["bam-files"].as<std::vector<path> >();
    }

    // Acquire path to genome file
    if (vm.count("genome-file")) {
        genomeFile = vm["genome-file"].as<path>();
    }

    // Test if provided genome exists
    if (!exists(genomeFile) && !symbolic_link_exists(genomeFile)) {
        BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
                    "Could not find genome file at: ") + genomeFile.string()));
    }

    // Glob the input bam files
    std::vector<path> transformedBams = Prepare::globFiles(bamFiles);

    auto_cpu_timer timer(1, "\nPortcullis completed.\nTotal runtime: %ws\n\n");        

    cout << "Running full portcullis pipeline" << endl
         << "--------------------------------" << endl << endl;


    if (!exists(outputDir)) {
        if (!create_directory(outputDir)) {
            BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
                    "Could not create output directory: ") + outputDir.string()));
        }
    }
    
    // ************ Prepare input data (BAMs + genome) ***********
    
    cout << "Preparing input data (BAMs + genome)" << endl
         << "----------------------------------" << endl << endl;
    
    path prepDir = path(outputDir.string() + "/prepare");

    // Create the prepare class
    Prepare prep(prepDir, strandednessFromString(strandSpecific), force, useLinks, useCsi, threads, verbose);
    
    // Prep the input to produce a usable indexed and sorted bam plus, indexed
    // genome and queryable coverage information
    prep.prepare(transformedBams, genomeFile);

    // Output any remaining details
    prep.outputDetails();
    
    
    // ************ Identify all junctions and calculate metrics ***********
    
    cout << "Identifying junctions and calculating metrics" << endl
         << "---------------------------------------------" << endl << endl;
    
    path juncDir = outputDir.string() + "/all_junctions";
    
    // Identify junctions and calculate metrics
    JunctionBuilder jb(prepDir.string(), juncDir.string(), "portcullis_all");
    jb.setThreads(threads);
    jb.setExtra(false);
    jb.setSeparate(false);
    jb.setVerbose(verbose);
    
    jb.process();
    

    // ************ Use default filtering strategy *************
    
    cout << "Filtering junctions" << endl
         << "-------------------" << endl << endl;
    
    path filtDir = outputDir.string() + "/filtered_junctions";
    path juncTab = juncDir.string() + "/portcullis_all.junctions.tab";
    
    JunctionFilter filter(juncTab, JunctionFilter::defaultFilterFile, filtDir, "portcullis_filtered");
    filter.setVerbose(verbose);
    filter.filter();

    
    // *********** BAM filter *********
    
    cout << "Filtering BAMs" << endl
         << "--------------" << endl << endl;
    
    path filtJuncTab = path(filtDir.string() + "/portcullis_filtered.junctions.tab");
    path bamFile = path(prepDir.string() + "/portcullis.sorted.alignments.bam");
    path filteredBam = path(outputDir.string() + "/portcullis.filtered.bam");
    
    BamFilter bamFilter(filtJuncTab.string(), bamFile.string(), filteredBam.string(), verbose);
    bamFilter.setStrandSpecific(strandednessFromString(strandSpecific));
    bamFilter.filter();

    return 0;        
}

/**
 * Start point for portcullis.
 */
int main(int argc, char *argv[]) {
    
    try {
        // Portcullis args
        string modeStr;
        std::vector<string> others;
        bool verbose;
        bool version;
        bool help;

        // Declare the supported options.
        po::options_description generic_options(helpHeader(), 100);
        generic_options.add_options()
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print extra information")
                ("version,V", po::bool_switch(&version)->default_value(false), "Print version string")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("mode", po::value<string>(&modeStr), "Portcullis mode.")
                ("others", po::value< std::vector<string> >(&others), "Other options.")
                ;

        // Positional options
        po::positional_options_description p;
        p.add("mode", 1);
        p.add("others", 100);
        
        // Combine non-positional options
        po::options_description cmdline_options;
        cmdline_options.add(generic_options).add(hidden_options);

        // Parse command line
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).allow_unregistered().run(), vm);
        po::notify(vm);

        // Output help information the exit if requested
        if (argc == 1 || argc == 2 && help) {
            cout << generic_options << endl;
            return 1;
        }

        // Always output version information but exit if version info was all user requested
        
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "Portcullis"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.11.X"
#endif
        
        // End if version was requested.
        if (version) {    
            cout << PACKAGE_NAME << " " << PACKAGE_VERSION << endl;
            return 0;
        }
        
        cout << "Portcullis V" << PACKAGE_VERSION << endl << endl;
        
        PortcullisFS fs(argv[0], PACKAGE_VERSION);
        
        // End if verbose was requested at this level, outputting file system details.
        if (verbose) {
            cout << endl 
                 << "Project filesystem" << endl 
                 << "------------------" << endl
                 << fs << endl;
            //return 0;
        }        
        
        // If we've got this far parse the command line properly
        Mode mode = parseMode(modeStr);
        
        const int modeArgC = argc-1;
        char** modeArgV = argv+1;
        
        // Set static variables in downstream subtools so they know where to get their resources from
        BamHelper::samtoolsExe = fs.getSamtoolsExe();
        JunctionFilter::defaultFilterFile = path(fs.getEtcDir().string() + "/default_filter.json");
        JunctionFilter::filterJuncsPy = fs.getFilterJuncsPy();
        JunctionSystem::version = fs.getVersion();
        
        if (mode == PREP) {
            Prepare::main(modeArgC, modeArgV);
        }
        else if(mode == JUNC) {
            JunctionBuilder::main(modeArgC, modeArgV);
        }
        else if (mode == FILTER) {
            JunctionFilter::main(modeArgC, modeArgV);
        }
        else if (mode == BAM_FILT) {
            BamFilter::main(modeArgC, modeArgV);
        }
        else if (mode == FULL) {
            mainFull(modeArgC, modeArgV);
        }
        else if (mode == CLUSTER) {
            Cluster::main(modeArgC, modeArgV);            
        }
        else if (mode == TRAIN) {
            Train::main(modeArgC, modeArgV);            
        }
        else if (mode == TEST) {
            Test::main(modeArgC, modeArgV);            
        }
        else {
            BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
                    "Unrecognised portcullis mode: ") + modeStr));
        }
                
    } catch(po::error& e) { 
        cerr << "Error: Parsing Command Line: " << e.what() << endl; 
        return 1; 
    } 
    catch (boost::exception &e) { 
        cerr << boost::diagnostic_information(e); 
        return 4;
    } catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 5;
    } catch (const char* msg) {
        cerr << "Error: " << msg << endl;
        return 6;
    } catch (...) {
        cerr << "Error: Exception of unknown type!" << endl;
        return 7;
    }

    return 0;
}

