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

#include "samtools_helper.hpp"
#include "intron.hpp"
#include "junction.hpp"
#include "junction_system.hpp"
#include "prepare.hpp"
using portcullis::BamAlignment;
using portcullis::BamReader;
using portcullis::GenomeMapper;
using portcullis::Intron;
using portcullis::Junction;
using portcullis::JunctionSystem;

#include "junction_builder.hpp"

portcullis::JunctionBuilder::JunctionBuilder(const path& _prepDir, const path& _outputDir, string _outputPrefix, uint16_t _threads, bool _fast, bool _verbose) {
    prepData = make_shared<PreparedFiles>(_prepDir);
    outputDir = _outputDir;
    outputPrefix = _outputPrefix;
    threads = _threads;
    fast = _fast;
    verbose = _verbose;        

    if (verbose) {
        cout << "Initialised Portcullis instance with settings:" << endl
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

    // Test if we have all the requried data
    if (!prepData->valid()) {
        BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Prepared data is not complete: ") + prepData->getPrepDir().string()));
    }

    if (verbose) {
        cout << "done." << endl << endl
             << "Loading settings stored in prep data ... ";
    }

    // Loading settings stored in prep data
    strandSpecific = prepData->loadSettings();

    if (verbose) {
        cout << "done." << endl
             << " - Strand specific input data: " << SSToString(strandSpecific) << endl << endl;
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

    const path sortedBamFile = prepData->getSortedBamFilePath();

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

    cout << endl << "Processing alignments and calculating metrics for reference sequence: " << endl;

    // The contents inside the pointer will automatically alter as reader.next() is called
    BamAlignment al;
    
    while(reader.next(al)) {

        while (junctionSystem.size() > 0 && lastCalculatedJunctionIndex < junctionSystem.size() && 
                al.getPosition() > junctionSystem.getJunctionAt(lastCalculatedJunctionIndex)->getIntron()->end) {

            JunctionPtr j = junctionSystem.getJunctionAt(lastCalculatedJunctionIndex);            

            /*size_t nbJABefore = j->getNbJunctionAlignmentFromVector();
            BamAlignmentPtr j0 = j->getFirstAlignment();
            size_t refBefore = j0.use_count() - 1;
            size_t alignSize = sizeof(*j0);
            size_t juncSizeBefore = sizeof(*j);*/

            j->calcMetrics();
            j->clearAlignments();
            lastCalculatedJunctionIndex++;

            /*size_t nbJAAfter = j->getNbJunctionAlignmentFromVector();
            size_t refAfter = j0.use_count() - 1;
            size_t juncSizeAfter = sizeof(*j);

            cout << "Junction " << lastCalculatedJunctionIndex << " shut down.  Intron start: " << j->getIntron()->start 
                    << ". Alignments: " << nbJABefore << "/" << nbJAAfter 
                    << ". Refcount: " << refBefore << "/" << refAfter 
                    << ". Size: " << juncSizeBefore << "/" << juncSizeAfter
                    << ". Alignment size: " << alignSize << "." << endl;*/

        }

        if (lastRefId == -1 || al.getReferenceId() != lastRefId) {
            cout << " - " << refs[al.getReferenceId()].name << endl;
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

        //cout << "After junction add (outer): " << bap.use_count() << endl;
    }

    while (junctionSystem.size() > 0 && lastCalculatedJunctionIndex < junctionSystem.size()) {

        JunctionPtr j = junctionSystem.getJunctionAt(lastCalculatedJunctionIndex);
        j->calcMetrics();
        j->clearAlignments();
        lastCalculatedJunctionIndex++;
    }

    reader.close();
    unsplicedWriter.close();
    splicedWriter.close();


    cout << endl << "Calculating metrics that require analysis of reference genome";
    if (verbose) {
        cout << endl;
    }
    else {
        cout << "...";
        cout.flush();
    }

    GenomeMapper gmap(prepData->getGenomeFilePath());
    gmap.loadFastaIndex();
    junctionSystem.scanReference(&gmap, refs, verbose);

    if (!verbose) {
        cout << " done" << endl;
    }

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
        junctionSystem.findFlankingAlignments(unsplicedFile, strandSpecific);
        
        cout << endl << "Calculating unspliced alignment coverage around junctions..." << endl;
        junctionSystem.calcCoverage(unsplicedFile, strandSpecific);
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
    string unsplicedIndexCmd = SamtoolsHelper::createIndexBamCmd(unsplicedFile);                

    int unsExitCode = system(unsplicedIndexCmd.c_str());                    

    cout << "done." << endl;
    cout << " - spliced alignments ... ";
    cout.flush();

    // Create BAM index
    string splicedIndexCmd = SamtoolsHelper::createIndexBamCmd(splicedFile);                

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


    auto_cpu_timer timer(1, "\nPortcullis junc completed.\nTotal runtime: %ws\n\n");        

    cout << "Running portcullis in junction builder mode" << endl
         << "------------------------------------------" << endl << endl;

    // Do the work ...
    JunctionBuilder jb(prepDir, outputDir, outputPrefix, threads, fast, verbose);
    jb.setSamtoolsExe(SamtoolsHelper::samtoolsExe);
    jb.process();

    return 0;
}

