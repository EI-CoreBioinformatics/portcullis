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

#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
using std::boolalpha;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using namespace boost::filesystem;
using boost::filesystem::path;
namespace po = boost::program_options;

#include "bam/bam_reader.hpp"
#include "bam/bam_writer.hpp"
#include "bam/depth_parser.hpp"
#include "bam/genome_mapper.hpp"
using namespace portcullis::bam;

#include "intron.hpp"
#include "junction.hpp"
#include "junction_system.hpp"
#include "prepare.hpp"
using portcullis::Intron;
using portcullis::Junction;
using portcullis::JunctionSystem;

#include "junction_builder.hpp"

portcullis::JunctionBuilder::JunctionBuilder(const path& _prepDir, const path& _outputDir, string _outputPrefix, uint16_t _threads, bool _fast, bool _verbose) {
    prepData = PreparedFiles(_prepDir);
    outputDir = _outputDir;
    outputPrefix = _outputPrefix;
    threads = _threads;
    fast = _fast;
    verbose = _verbose;        

    if (verbose) {
        cout << "Initialised Portcullis instance with settings:" << endl
             << " - Prep data dir: " << prepData.getPrepDir() << endl
             << " - Output directory: " << outputDir << endl
             << " - Output file name prefix: " << outputPrefix << endl
             << " - Threads: " << threads << endl 
             << " - Fast mode: " << boolalpha << fast << endl << endl;            
    }

    if (verbose) {
        cout << "Ensuring output directory exists ... ";
        cout.flush();
    }

    if (!bfs::exists(outputDir)) {
        if (!bfs::create_directory(outputDir)) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not create output directory: ") + outputDir.string()));
        }
    }

    if (verbose) {
        cout << "done." << endl << endl
             << "Checking prepared data ... ";
        cout.flush();
    }
    
    // Loading settings stored in prep data
    settings = prepData.loadSettings();

    // Test if we have all the required data
    if (!prepData.valid(settings.useCsi)) {
        BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Prepared data is not complete: ") + prepData.getPrepDir().string()));
    }

    if (verbose) {
        cout << "done." << endl
             << " - Strand specific input data: " << SSToString(settings.ss) << endl
             << " - Indexing format: " << (settings.useCsi ? "CSI" : "BAI") << endl << endl;
    }
}
    
portcullis::JunctionBuilder::~JunctionBuilder() {

    splicedAlignmentMap.clear();
}
    
    
    
/**
 * Populates the set of distinct junctions.  
 * 
 * Also outputs all the unspliced alignments to a separate file if requested
 */
void portcullis::JunctionBuilder::process() {

    auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");

    const path sortedBamFile = prepData.getSortedBamFilePath();

    BamReader reader(sortedBamFile, 1);
    reader.open();

    vector<RefSeq> refs = reader.getRefs();
    junctionSystem.setRefs(refs);

    cout << " - Reading alignments from: " << sortedBamFile << endl;

    path unsplicedFile = getUnsplicedBamFile();
    path splicedFile = getSplicedBamFile();

    BamWriter unsplicedWriter(unsplicedFile);
    BamWriter splicedWriter(splicedFile);

    cout << " - Saving unspliced alignments to: " << unsplicedFile << endl;
    unsplicedWriter.open(reader.getHeader());

    cout << " - Saving spliced alignments to: " << splicedFile << endl;
    splicedWriter.open(reader.getHeader());

    uint64_t splicedCount = 0;
    uint64_t unsplicedCount = 0;
    uint64_t sumQueryLengths = 0;
    int32_t minQueryLength = 100000;
    int32_t maxQueryLength = 0;

    int32_t lastRefId = -1;
    int32_t lastCalculatedJunctionIndex = 0;
    int32_t lastSeqCount = 0;
    int32_t chunkSize = 0;
    int32_t nextTarget = 0;
    int32_t percentComplete = 0;

    cout << endl << "Loading reference index ... ";
    cout.flush();
    
    GenomeMapper gmap(prepData.getGenomeFilePath());
    gmap.loadFastaIndex();
    
    cout << "done." << endl << endl;
        
    cout << "Processing alignments and calculating metrics for reference sequence: " << endl;

    // The contents inside the pointer will automatically alter as reader.next() is called
    
    while(reader.next()) {

        const BamAlignment& al = reader.current();
        
        while (junctionSystem.size() > 0 && lastCalculatedJunctionIndex < junctionSystem.size() && 
                (al.getReferenceId() != lastRefId || al.getPosition() > junctionSystem.getJunctionAt(lastCalculatedJunctionIndex)->getIntron()->end)) {

            JunctionPtr j = junctionSystem.getJunctionAt(lastCalculatedJunctionIndex);            

            //cout << "Processing junction: " << lastCalculatedJunctionIndex << " - " << j->getIntron()->toString() << endl;
            j->calcMetrics();
            j->processJunctionWindow(gmap);
            j->clearAlignments();
            
            lastCalculatedJunctionIndex++;
            lastSeqCount++;
        }

        if (lastRefId == -1 || al.getReferenceId() != lastRefId) {
            
            // Just check we haven't got any strange alignments that are not associated with ref seqs.
            // End if we do
            if (al.getReferenceId() >= refs.size()) {
                break;
            }
            
            if (lastRefId > -1) {
                cout << "100%\t" << lastSeqCount << " potential junctions found." << endl;
                nextTarget = 0;
                percentComplete = 0;
            }
            
            cout << " - " << refs[al.getReferenceId()].name << "\t0% ... ";
            cout.flush();
            
            chunkSize = refs[al.getReferenceId()].length / 10;
            nextTarget += chunkSize;
            percentComplete += 10;
            
            lastSeqCount = 0;
        }

        lastRefId = al.getReferenceId();

        int32_t len = al.getLength();
        minQueryLength = min(minQueryLength, len);
        maxQueryLength = max(maxQueryLength, len);

        sumQueryLengths += len;

        //cout << "Before junction add (outer): " << bap.use_count() << endl;            
        if (junctionSystem.addJunctions(al, false)) {//strandSpecific)) {
            splicedWriter.write(al);
            splicedCount++;

            // Record alignment name in map
            size_t code = std::hash<string>()(al.deriveName());
            splicedAlignmentMap[code]++;
        }
        else {
            unsplicedWriter.write(al);
            unsplicedCount++;
        }

        while (al.getPosition() > nextTarget) {
            cout << percentComplete << "% ... ";
            cout.flush();
            
            nextTarget += chunkSize;
            percentComplete += 10; 
        }
    }

    while (junctionSystem.size() > 0 && lastCalculatedJunctionIndex < junctionSystem.size()) {

        JunctionPtr j = junctionSystem.getJunctionAt(lastCalculatedJunctionIndex);
        j->calcMetrics();
        j->processJunctionWindow(gmap);
        j->clearAlignments();
        lastCalculatedJunctionIndex++;
        lastSeqCount++;
    }
    
    cout << "100%\t" << lastSeqCount << " potential junctions found." << endl;

    reader.close();
    unsplicedWriter.close();
    splicedWriter.close();


    cout << endl << "Calculating junctions stats that require comparisons with other junctions...";
    cout.flush();
    junctionSystem.calcJunctionStats();
    junctionSystem.calcMultipleMappingStats(splicedAlignmentMap);
    cout << " done" << endl;

    if (fast) {
        //cout << "skipped due to user request to run in fast mode" << endl << endl;
    }
    else {
        // Count the number of alignments found in upstream and downstream flanking 
        // regions for each junction
        cout << endl << endl << "Analysing unspliced alignments around junctions:" << endl;
        junctionSystem.findFlankingAlignments(unsplicedFile, settings.ss);
        
        cout << endl << "Calculating unspliced alignment coverage around junctions..." << endl;
        junctionSystem.calcCoverage(unsplicedFile, settings.ss);
    }

    // Calculate some stats
    uint64_t totalAlignments = splicedCount + unsplicedCount;
    double meanQueryLength = (double)sumQueryLengths / (double)totalAlignments;
    junctionSystem.setQueryLengthStats(minQueryLength, meanQueryLength, maxQueryLength);

    cout << "Processed " << totalAlignments << " alignments." << endl
         << " - Alignment query length statistics: min: " << minQueryLength << "; mean: " << meanQueryLength << "; max: " << maxQueryLength << ";" << endl
         << " - Found " << junctionSystem.size() << " junctions from " << splicedCount << " spliced alignments." << endl
         << " - Found " << unsplicedCount << " unspliced alignments." << endl << endl;

    cout << "Indexing:" << endl;
    cout << " - unspliced alignments ... ";

    // Create BAM index
    string unsplicedIndexCmd = BamHelper::createIndexBamCmd(unsplicedFile, settings.useCsi);                

    int unsExitCode = system(unsplicedIndexCmd.c_str());                    

    cout << "done." << endl;
    cout << " - spliced alignments ... ";
    cout.flush();

    // Create BAM index
    string splicedIndexCmd = BamHelper::createIndexBamCmd(splicedFile, settings.useCsi);                

    int sExitCode = system(splicedIndexCmd.c_str());                    

    cout << "done." << endl;

    cout << endl << "Saving junctions: " << endl;
    junctionSystem.saveAll(path(outputDir.string() + "/" + outputPrefix));
}
   
int portcullis::JunctionBuilder::main(int argc, char *argv[]) {

    // Portcullis args
    string prepDir;
    string outputDir;
    string outputPrefix;
    uint16_t threads;
    bool fast;
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(helpMessage(), 100);
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


    auto_cpu_timer timer(1, "\nPortcullis junc completed.\nTotal runtime: %ws\n\n");        

    cout << "Running portcullis in junction builder mode" << endl
         << "------------------------------------------" << endl << endl;

    // Do the work ...
    JunctionBuilder jb(prepDir, outputDir, outputPrefix, threads, fast, verbose);
    jb.setSamtoolsExe(BamHelper::samtoolsExe);
    jb.process();

    return 0;
}

