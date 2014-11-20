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
#include <iostream>
#include <vector>
using std::boolalpha;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
namespace po = boost::program_options;

#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;

#include "genome_mapper.hpp"
#include "intron.hpp"
#include "junction.hpp"
#include "junction_system.hpp"
#include "prepare.hpp"
using portculis::GenomeMapper;
using portculis::Intron;
using portculis::Junction;
using portculis::JunctionSystem;



namespace portculis {

const string DEFAULT_JUNC_OUTPUT_DIR = "portculis_junc_out";
const string DEFAULT_JUNC_OUTPUT_PREFIX = "portculis";
const uint16_t DEFAULT_JUNC_THREADS = 1;

typedef boost::error_info<struct JunctionBuilderError,string> JunctionBuilderErrorInfo;
struct JunctionBuilderException: virtual boost::exception, virtual std::exception { };

class JunctionBuilder {
private:

    // Can set these from the outside via the constructor
    PreparedFiles* prepData;
    string outputDir;
    string outputPrefix;
    bool strandSpecific;
    uint16_t threads;
    bool fast;
    bool verbose;
    
    // Sam header and refs info from the input bam
    SamHeader header;
    RefVector refs;
    
    // The set of distinct junctions found in the BAM file
    JunctionSystem junctionSystem;
    

protected:
    
    string getUnsplicedBamFile() {
        return outputDir + "/" + outputPrefix + ".unspliced.bam";
    }
    
    string getAssociatedIndexFile(string bamFile) {
        return string(bamFile) + ".bti";
    }
    
    /**
     * Populates the set of distinct junctions.  
     * 
     * Also outputs all the unspliced alignments to a separate file if requested
     */
    void separateSplicedAlignments() {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
        
        BamReader reader;
        
        const string sortedBamFile = prepData->getSortedBamFilePath();
        
        if (!reader.Open(sortedBamFile)) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not open BAM reader for input: ") + sortedBamFile));
        }
        
        // Sam header and refs info from the input bam
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

        junctionSystem.setRefs(refs);
       
        cout << " - Reading alignments from: " << sortedBamFile << endl;
        
        string indexFile = prepData->getBamIndexFilePath();
        
        // Opens the index for this BAM file
        if ( !reader.OpenIndex(indexFile) ) {            
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not open index for BAM: ") + indexFile));             
        }
        
        cout << " - Using BAM index: " << indexFile << endl;
        
        BamWriter unsplicedWriter;
        string unsplicedFile = getUnsplicedBamFile();

        if (!unsplicedWriter.Open(unsplicedFile, header, refs)) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not open BAM writer for non-spliced file: ") + unsplicedFile));
        }

        cout << " - Saving unspliced alignments to: " << unsplicedFile << endl;
        
        BamAlignment al;
        uint64_t splicedCount = 0;
        uint64_t unsplicedCount = 0;
        uint64_t sumQueryLengths = 0;
        int32_t minQueryLength = 100000;
        int32_t maxQueryLength = 0;
        cout << " - Processing alignments ... ";
        cout.flush();
        while(reader.GetNextAlignment(al))
        {
            int32_t len = al.Length;
            minQueryLength = min(minQueryLength, len);
            maxQueryLength = max(maxQueryLength, len);
            
            sumQueryLengths += len;
            
            if (junctionSystem.addJunctions(al, false)) {//strandSpecific)) {
                splicedCount++;
            }
            else {
                unsplicedWriter.SaveAlignment(al);
                unsplicedCount++;
            }
        }
        
        reader.Close();
        unsplicedWriter.Close();
        cout << "done." << endl;
        
        // Calculate some stats
        uint64_t totalAlignments = splicedCount + unsplicedCount;
        double meanQueryLength = (double)sumQueryLengths / (double)totalAlignments;
        junctionSystem.setQueryLengthStats(minQueryLength, meanQueryLength, maxQueryLength);
        
        cout << " - Processed " << totalAlignments << " alignments." << endl
             << " - Alignment query length statistics: min: " << minQueryLength << "; mean: " << meanQueryLength << "; max: " << maxQueryLength << ";" << endl
             << " - Found " << junctionSystem.size() << " junctions from " << splicedCount << " spliced alignments." << endl
             << " - Found " << unsplicedCount << " unspliced alignments." << endl;
        
        
        BamReader indexReader;
        if (!reader.Open(unsplicedFile)) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not open bam reader for unspliced alignments file: ") + unsplicedFile));
        }
        // Sam header and refs info from the input bam
        SamHeader header = reader.GetHeader();
        RefVector refs = reader.GetReferenceData();

        // Opens the index for this BAM file
        string unsplicedIndexFile = getAssociatedIndexFile(unsplicedFile);
        if ( !reader.OpenIndex(unsplicedIndexFile) ) {            
            if ( !reader.CreateIndex(BamIndex::BAMTOOLS) ) {
                BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                        "Error creating BAM index for unspliced alignments file: ") + unsplicedIndexFile));
            }            
        }
    }
    

public:

    
    JunctionBuilder(string _prepDir, string _outputDir, string _outputPrefix, uint16_t _threads, bool _fast, bool _verbose) {
        prepData = new PreparedFiles(_prepDir);
        outputDir = _outputDir;
        outputPrefix = _outputPrefix;
        threads = _threads;
        fast = _fast;
        verbose = _verbose;        
        
        if (verbose) {
            cout << "Initialised Portculis instance with settings:" << endl
                 << " - Prep data dir: " << prepData->getPrepDir() << endl
                 << " - Output directory: " << outputDir << endl
                 << " - Output file name prefix: " << outputPrefix << endl
                 << " - Threads: " << threads << endl 
                 << " - Fast mode: " << boolalpha << fast << endl << endl;            
        }
        
        if (verbose) {
            cout << "Ensuring output directory exists ... ";
            cout.flush();
        }
        
        if (!exists(outputDir)) {
            if (!create_directory(outputDir)) {
                BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                        "Could not create output directory: ") + outputDir));
            }
        }
        
        if (verbose) {
            cout << "done." << endl << endl
                 << "Checking prepared data ... ";
            cout.flush();
        }
        
        // Test if we have all the requried data
        if (!prepData->valid()) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                        "Prepared data is not complete: ") + prepData->getPrepDir()));
        }
        
        if (verbose) {
            cout << "done." << endl << endl
                 << "Loading settings stored in prep data ... ";
        }
        
        // Loading settings stored in prep data
        strandSpecific = prepData->loadSettings();
        
        if (verbose) {
            cout << "done." << endl
                 << " - Strand specific input data: " << boolalpha << strandSpecific << endl << endl;
        }
    }
    
    virtual ~JunctionBuilder() {

        if (prepData != nullptr) {
            delete prepData;
        }        
    }
    

    void process() {
       
        // Collect junctions from BAM file (also outputs unspliced alignments
        // to a separate file)
        cout << "Stage 1: Separating spliced alignments:" << endl;
        separateSplicedAlignments();        
        
        // Acquires donor / acceptor info from indexed genome file
        cout << "Stage 2: Scanning reference sequences:" << endl;
        GenomeMapper gmap(prepData->getGenomeFilePath());
        gmap.loadFastaIndex();
        junctionSystem.scanReference(&gmap, refs);
        
        if (fast) {
            cout << "Stage 3: skipped due to user request to run in fast mode" << endl << endl;
        }
        else {
            // Count the number of alignments found in upstream and downstream flanking 
            // regions for each junction
            cout << "Stage 3: Analysing unspliced alignments around junctions:" << endl;
            junctionSystem.findFlankingAlignments(getUnsplicedBamFile(), strandSpecific);
        }

        cout << "Stage 4: Calculating unspliced alignment coverage around junctions:" << endl;
        junctionSystem.calcCoverage(getUnsplicedBamFile(), strandSpecific);
            
        cout << "Stage 5: Calculating junction status flags:" << endl;
        junctionSystem.calcJunctionStats();
        
        // Calculate all remaining metrics
        cout << "Stage 6: Calculating remaining junction metrics:" << endl;
        junctionSystem.calcAllRemainingMetrics();
        
        cout << "Stage 7: Outputting junction information:" << endl;
        junctionSystem.saveAll(outputDir + "/" + outputPrefix);
    }
    
    static string helpMessage() {
        return string("\nPortculis Junction Builder Mode Help.\n\n") +
                      "Usage: portculis junc [options] <prep_data_dir> \n\n" +
                      "Run \"portculis prep ...\" to generate data suitable for junction finding before running \"portculis junc ...\"\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portculis args
        string prepDir;
        string outputDir;
        string outputPrefix;
        uint16_t threads;
        bool fast;
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output_dir,o", po::value<string>(&outputDir)->default_value(DEFAULT_JUNC_OUTPUT_DIR), 
                    "Output directory for files generated by this program.")
                ("output_prefix,p", po::value<string>(&outputPrefix)->default_value(DEFAULT_JUNC_OUTPUT_PREFIX), 
                    "File name prefix for files generated by this program.")
                ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_JUNC_THREADS),
                    (string("The number of threads to use.  Default: ") + lexical_cast<string>(DEFAULT_JUNC_THREADS)).c_str())
                ("fast,f", po::bool_switch(&fast)->default_value(false),
                    "If running in fast mode then metrics 10 (up), 11 (down) are not calculated")
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                    "Print extra information")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("prep_data_dir,i", po::value<string>(&prepDir), "Path to directory containing prepared data.")
                ;

        // Positional option for the input bam file
        po::positional_options_description p;
        p.add("prep_data_dir", 1);
        
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
        
        // Acquire path to genome file
        if (vm.count("prep_data_dir")) {
            prepDir = vm["prep_data_dir"].as<string>();
        }
        
        
        auto_cpu_timer timer(1, "\nPortculis junc completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portculis in junction builder mode" << endl
             << "------------------------------------------" << endl << endl;
        
        // Do the work ...
        JunctionBuilder(prepDir, outputDir, outputPrefix, threads, fast, verbose).process();
        
        return 0;
    }
};
}
