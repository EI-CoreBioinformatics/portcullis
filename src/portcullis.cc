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

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "junction_builder.hpp"
#include "genome_mapper.hpp"
#include "prepare.hpp"
#include "filter.hpp"
#include "cluster.hpp"
#include "train.hpp"
#include "test.hpp"
using portcullis::JunctionBuilder;
using portcullis::Prepare;
using portcullis::Filter;
using portcullis::Cluster;
using portcullis::Train;
using portcullis::Test;

#include "train.hpp"

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
                  " - full    - Runs prep, junc, filter as a complete pipeline\n" +
                  " - cluster - Clusters potential junctions to help distinguish real junctions for false\n\n" +
                  "\nAvailable options";
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
        po::options_description generic_options(helpHeader());
        generic_options.add_options()
                ("version", po::bool_switch(&version)->default_value(false), "Print version string")
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

        // Output version information then exit if requested
        if (version) {
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "Portcullis"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.2.1"
#endif
            cout << PACKAGE_NAME << " V" << PACKAGE_VERSION << endl;
            return 0;
        }
        
        
        Mode mode = parseMode(modeStr);
        
        const int modeArgC = argc-1;
        char** modeArgV = argv+1;
        
        if (mode == PREP) {
            Prepare::main(modeArgC, modeArgV);
        }
        else if(mode == JUNC) {
            JunctionBuilder::main(modeArgC, modeArgV);
        }
        else if (mode == FILTER) {
            Filter::main(modeArgC, modeArgV);
        }
        else if (mode == FULL) {
            
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
                
    } catch (boost::exception &e) { 
        std::cerr << boost::diagnostic_information(e); 
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

