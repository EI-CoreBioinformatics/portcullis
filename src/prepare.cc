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

#include <sys/ioctl.h>
#include <climits>
#include <glob.h>
#include <fstream>
#include <string>
#include <memory>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::shared_ptr;
using std::make_shared;
using std::vector;
using std::istream_iterator;

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
#include <portcullis/bam/genome_mapper.hpp>
#include <portcullis/portcullis_fs.hpp>
using portcullis::PortcullisFS;
using namespace portcullis::bam;

#include "prepare.hpp"

bool portcullis::PreparedFiles::valid(bool useCsi) const {
        
    if (!bfs::exists(getSortedBamFilePath()) && !bfs::symbolic_link_exists(getSortedBamFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find sorted BAM files at: ") + getSortedBamFilePath().string()));
    }

    if (!bfs::exists(getBamIndexFilePath(useCsi)) && !bfs::symbolic_link_exists(getBamIndexFilePath(useCsi))) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find BAM index at: ") + getBamIndexFilePath(useCsi).string()));
    }

    if (!bfs::exists(getGenomeFilePath()) && !bfs::symbolic_link_exists(getGenomeFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find genome file at: ") + getGenomeFilePath().string()));
    }

    if (!bfs::exists(getGenomeIndexFilePath()) && !bfs::symbolic_link_exists(getGenomeIndexFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find genome index at: ") + getGenomeIndexFilePath().string()));
    }

    return true;
}
    
void portcullis::PreparedFiles::clean() {

    bfs::remove(getUnsortedBamFilePath());
    bfs::remove(getSortedBamFilePath());
    bfs::remove(getBamIndexFilePath(false));
    bfs::remove(getBamIndexFilePath(true));
    bfs::remove(getGenomeFilePath());
    bfs::remove(getGenomeIndexFilePath());
    bfs::remove(getBcfFilePath());
    bfs::remove(getBcfIndexFilePath());
}
    
    
portcullis::Prepare::Prepare(const path& _outputPrefix, Strandedness _strandSpecific, bool _force, bool _useLinks, bool _useCsi, uint16_t _threads, bool _verbose) {
    output = make_shared<PreparedFiles>(_outputPrefix);
    strandSpecific = _strandSpecific;
    force = _force;
    useLinks = _useLinks;
    useCsi = _useCsi;
    threads = _threads;
    verbose = _verbose;

    if (verbose) {
        cout << "Configured portcullis prep to use the following settings: " << endl
             << " - Output directory: " << output->getPrepDir() << endl
             << " - Strand specific library: " << strandednessToString(strandSpecific) << endl
             << " - Force prep (cleans output directory): " << boolalpha << force << endl
             << " - Use symbolic links instead of copy where possible: " << boolalpha << useLinks << endl
             << " - Indexing type: " << (useCsi ? "CSI" : "BAI") << endl
             << " - Threads (for sorting BAM): " << threads << endl << endl;
    }

    if (force) {
        cout << "Cleaning output dir " << _outputPrefix << " ... ";
        cout.flush();
        output->clean();
        cout << "done." << endl << endl;
    }
}
    
bool portcullis::Prepare::copy(const path& from, const path& to, const string& msg, const bool requireInputFileExists) {

    
    // Check if output file already exists (if so don't do anything)
    if (bfs::exists(to) || bfs::symbolic_link_exists(to)) {
        cout << "Prepped " << msg << " file detected: " << to << endl;            
    }
    // Check if input file is required to exist (if so either symlink or copy)
    else if (requireInputFileExists || bfs::exists(from) || bfs::symbolic_link_exists(from)) {

        if (useLinks) {
            bfs::create_symlink(bfs::canonical(from), to);
            cout << "Created symlink from " << from << " to " << to << endl;
        }
        else {
        
            auto_cpu_timer timer(1, string(" - Copy ") + msg + " - Wall time taken: %ws\n\n");
            
            cout << "Copying from " << from << " to " << to << " ... ";
            cout.flush();
            
            std::ifstream  src(from.string(), std::ios::binary);
            std::ofstream  dst(to.string(),   std::ios::binary);

            dst << src.rdbuf();

            src.close();
            dst.close();

            cout << "done." << endl;
        }
    }
    // If both input file and output file do not exist and it's not required to exists
    else if (!requireInputFileExists) {
        cout << "Existing " << msg << " not found.  Will create later." << endl;        
    }
    
    return bfs::exists(to) || bfs::symbolic_link_exists(to);
}
    
bool portcullis::Prepare::genomeIndex() {

    const path genomeFile = output->getGenomeFilePath();
    const path indexFile = output->getGenomeIndexFilePath();

    bool indexExists = bfs::exists(indexFile);

    if (indexExists) {
        cout << "Pre-indexed genome detected: " << indexFile << endl;            
    }
    else {

        auto_cpu_timer timer(1, " - Genome Index - Wall time taken: %ws\n\n");

        cout << "Indexing genome " << genomeFile << " ... ";
        cout.flush();

        // Create the index
        GenomeMapper(genomeFile).buildFastaIndex();

        cout << "done." << endl
             << "Genome index file created at: " << output->getGenomeIndexFilePath() << endl;                
    }

    return bfs::exists(indexFile);
}
    
    
/**
 * Merge together a set of BAM files, use the output prefix to construct a
 * file name for the merged file
 * @param bamFiles The set of bam files to merge
 * @return 
 */
bool portcullis::Prepare::bamMerge(vector<path> bamFiles) {

    
    path mergedBam = output->getSortedBamFilePath();        

    bool mergedBamExists = bfs::exists(mergedBam) || bfs::symbolic_link_exists(mergedBam);

    if (mergedBamExists) {
        cout << "Pre-merged BAM detected: " << mergedBam << endl;
    }
    else {

        auto_cpu_timer timer(1, " - BAM Merge - Wall time taken: %ws\n\n");

        cout << "Found " << bamFiles.size() << " BAM files." << endl;
        
        vector<path> mergeIn;
        
        // Sort the individual inputs if necessary
        uint16_t inCount = 1;
        for(auto& f : bamFiles) {
            path tempSorted = path(output->getPrepDir());
            tempSorted /= "temp" + lexical_cast<string>(inCount++) + ".bam";
            // Sort the data (if required, will auto-detect if necessary)
            if (!bamSort(f, tempSorted)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not sort: ") + output->getUnsortedBamFilePath().string()));
            }
            mergeIn.push_back(tempSorted);            
        }

        
        string mergeCmd = BamHelper::createMergeBamCmd(mergeIn, mergedBam, threads);

        cout << "Merging BAM using command \"" << mergeCmd << "\" ... ";
        cout.flush();

        int exitCode = system(mergeCmd.c_str());

        if (exitCode != 0 || !bfs::exists(mergedBam)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Failed to successfully merge: ") + mergedBam.string()));
        }

        cout << "done." << endl
             << "Merged BAM file created at: " << mergedBam << endl;
        
        cout << "Deleteing temporary unmerged files ...";
        cout.flush();
        for(auto& f : mergeIn) {
            bfs::remove(f);
        }
        cout << "done." << endl;
    }

    // Return true if the merged BAM exists now, which is should do
    return bfs::exists(mergedBam) || bfs::symbolic_link_exists(mergedBam);        
}

/**
 * Sorts the unsorted bam file if required or forced
 * @param inputBam
 * @return 
 */
bool portcullis::Prepare::bamSort(const path& input, const path& output) {

    const path unsortedBam = input;
    const path sortedBam = output;

    bool sortedBamExists = bfs::exists(sortedBam) || bfs::symbolic_link_exists(sortedBam);

    if (sortedBamExists) {            
        cout << "Prepped sorted BAM detected: " << sortedBam << endl;
    }
    else {

        if (BamHelper::isCoordSortedBam(unsortedBam.string()) && !force) {

            cout << "Provided BAM appears to be sorted already, just creating symlink instead." << endl;
            bfs::create_symlink(bfs::canonical(unsortedBam), sortedBam);            
            cout << "Created symlink from " << bfs::canonical(unsortedBam) << " to " << sortedBam << endl;
        }
        else {

            auto_cpu_timer timer(1, " - BAM Sort - Wall time taken: %ws\n\n");

            // Sort the BAM file by coordinate
            string sortCmd = BamHelper::createSortBamCmd(unsortedBam, sortedBam, false, threads, "1G");

            cout << "Sorting BAM using command \"" << sortCmd << "\" ... ";
            cout.flush();

            int exitCode = system(sortCmd.c_str());

            path badNameMergeFile = path(sortedBam.string() + ".bam");

            if (bfs::exists(badNameMergeFile) || bfs::symbolic_link_exists(badNameMergeFile)) {
                boost::filesystem::rename(badNameMergeFile, sortedBam);
            }

            if (exitCode != 0 || !bfs::exists(sortedBam) || !BamHelper::isCoordSortedBam(sortedBam)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Failed to successfully sort: ") + unsortedBam.string()));
            }

            cout << "done." << endl
                 << "Sorted BAM file created at: " << sortedBam << endl;            
        }
    }

    // Return true if the sorted BAM exists now, which is should do
    return bfs::exists(sortedBam) || bfs::symbolic_link_exists(sortedBam);
}
    
bool portcullis::Prepare::bamIndex(const bool copied) {

    const path sortedBam = output->getSortedBamFilePath();
    const path indexedFile = output->getBamIndexFilePath(useCsi);

    bool indexedBamExists = bfs::exists(indexedFile) || bfs::symbolic_link_exists(indexedFile);

    if (indexedBamExists && !copied) {
        if (verbose) cout << "Prepped indexed BAM detected: " << indexedFile << endl;
    }
    else if (!indexedBamExists) {

        auto_cpu_timer timer(1, " - BAM Index - Wall time taken: %ws\n\n");

        // Create BAM index
        string indexCmd = BamHelper::createIndexBamCmd(sortedBam, useCsi);                

        cout << "Indexing BAM using command \"" << indexCmd << "\" ... ";
        cout.flush();
        
        int exitCode = system(indexCmd.c_str());                    

        if (exitCode != 0 || !exists(output->getBamIndexFilePath(useCsi))) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Failed to successfully index: ") + sortedBam.string()));
        }

        cout << "done." << endl
             << "BAM index created at: " << output->getBamIndexFilePath(useCsi) << endl;
    }

    return bfs::exists(indexedFile) || bfs::symbolic_link_exists(indexedFile);
}
    
void portcullis::Prepare::prepare(vector<path> bamFiles, const path& originalGenomeFile) {

    // Copy / Symlink the genome file to the output dir
    if (!copy(originalGenomeFile, output->getGenomeFilePath(), "genome", true)) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not copy/symlink genome file to: ") + output->getGenomeFilePath().string()));
    }
    
    // Copy / Symlink the genome index file to the output dir, or create it
    if (!copy(originalGenomeFile.string() + FASTA_INDEX_EXTENSION, output->getGenomeFilePath().string() + FASTA_INDEX_EXTENSION, "genome index", false)) {
        if (!genomeIndex()) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not create genome index")));
        }
    }
    
    bool validIndexingMode = checkIndexMode(output->getGenomeFilePath(), useCsi);
    
    if (!validIndexingMode) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "User requested ") + (useCsi ? "CSI" : "BAI") + " indexing mode, however, genome file contains sequences too long to properly index using this method.  To continue, restart using the --use_csi option."));
    }
    
    const bool doMerge = bamFiles.size() > 1;        

    if (bamFiles.empty()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "No BAM files to process")));
    }

    path mergedBamFile = doMerge ? output->getUnsortedBamFilePath() : bamFiles[0];

    bool indexCopied = false;
    // Merge the bams to output a sorted bam file if required, otherwise just
    // copy / symlink the file provided
    if (doMerge) {
        if (!bamMerge(bamFiles)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not merge BAM files")));
        }
    }
    else {

        // Copy / Symlink the file to the output dir
        if (!copy(bamFiles[0], output->getUnsortedBamFilePath(), "BAM", true)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Could not copy/symlink BAM file to: ") + output->getUnsortedBamFilePath().string()));
        }
        
        // Copy / Symlink the index file to the output dir if it exists... otherwise we'll create it later
        indexCopied = copy(bamFiles[0].string() + (useCsi ? CSI_EXTENSION : BAI_EXTENSION), output->getBamIndexFilePath(useCsi), "BAM index", false);
        
        // Sort the data (if required, will auto-detect if necessary)
        if (!bamSort(output->getUnsortedBamFilePath(), output->getSortedBamFilePath())) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not sort: ") + output->getUnsortedBamFilePath().string()));
        }
        
        // Save disk space by deleting the unsorted BAM if present
        if (!bfs::symbolic_link_exists(output->getUnsortedBamFilePath()) && bfs::exists(output->getUnsortedBamFilePath())) {
            bfs::remove(output->getUnsortedBamFilePath());
        }
    }
    
    // Index the sorted file (if required)
    if (!bamIndex(indexCopied)) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Failed to index: ") + output->getSortedBamFilePath().string()));
    }
}

vector<path> portcullis::Prepare::globFiles(vector<path> input) {

    glob_t globbuf;

    // Translate glob patterns into real paths
    int i = 0;
    for(path g : input) {           
        glob(g.c_str(), i > 0 ? GLOB_TILDE | GLOB_APPEND : GLOB_TILDE, NULL, &globbuf);
        i++;
    }

    vector<path> transformedBams;
    for( int i = 0; i < globbuf.gl_pathc; ++i )
        transformedBams.push_back( path(globbuf.gl_pathv[i]) );

    if( globbuf.gl_pathc > 0 )
        globfree( &globbuf );

    return transformedBams;
}

bool portcullis::Prepare::checkIndexMode(const path& genomeIndexFile, const bool useCsi) {
    
    ifstream genomeIndex(genomeIndexFile.string());
    
    vector<string> myLines;
    std::copy(istream_iterator<string>(genomeIndex),
          istream_iterator<string>(),
          back_inserter(myLines));
    
    for (auto& l : myLines) {
        vector<string> parts;
        boost::split( parts, l, boost::is_any_of("\t"), boost::token_compress_on );
        
        if (parts.size() >= 2) {
            if (boost::lexical_cast<int>(parts[1]) >= LONG_MAX) {
                return false;
            }
        }    
    }
    
    return true;
}
    
int portcullis::Prepare::main(int argc, char *argv[]) {

    // Portcullis args
    vector<path> bamFiles;
    path genomeFile;
    path outputDir;
    bool force;
    string strandSpecific;
    bool copy;
    bool useCsi;
    uint16_t threads;
    bool verbose;
    bool help;
    
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    // Declare the supported options.
    po::options_description generic_options(helpMessage(), w.ws_col, (unsigned)((double)w.ws_col/1.7));
    generic_options.add_options()
            ("output,o", po::value<path>(&outputDir)->default_value(DEFAULT_PREP_OUTPUT_DIR), 
                (string("Output directory for prepared files. Default: ") + DEFAULT_PREP_OUTPUT_DIR).c_str())
            ("force", po::bool_switch(&force)->default_value(false), 
                "Whether or not to clean the output directory before processing, thereby forcing full preparation of the genome and bam files.  By default portcullis will only do what it thinks it needs to.")
            ("strandedness", po::value<string>(&strandSpecific)->default_value(strandednessToString(Strandedness::UNKNOWN)), 
                "Whether BAM alignments were generated using a strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads).  By default we assume the user does not know the strand specific protocol used for this BAM file.  This has the affect that strand information is derived from splice site information alone, assuming junctions are either canonical or semi-canonical in form.  Default: \"unknown\"")
            ("copy", po::bool_switch(&copy)->default_value(false), 
                "Whether to copy files from input data to prepared data where possible, otherwise will use symlinks.  Will require more time and disk space to prepare input but is potentially more robust.")
            ("use_csi,c", po::bool_switch(&useCsi)->default_value(false), 
                "Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it supports very long target sequences (probably not an issue unless you are working on huge genomes).  BAI has the advantage that it is more widely supported (useful for viewing in genome browsers).")
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
            ("bam-files", po::value< vector<path> >(&bamFiles), "Path to the BAM files to process.")
            ("genome-file", po::value<path>(&genomeFile), "Path to the genome file to process.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("genome-file", 1);
    p.add("bam-files", -1);

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

    // Acquire path to bam files
    if (vm.count("bam-files")) {
        bamFiles = vm["bam-files"].as<vector<path> >();
    }

    // Acquire path to genome file
    if (vm.count("genome-file")) {
        genomeFile = vm["genome-file"].as<path>();
    }

    // Test if provided genome exists
    if (!bfs::exists(genomeFile) && !bfs::symbolic_link_exists(genomeFile)) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not find genome file at: ") + genomeFile.string()));
    }

    if (bamFiles.empty()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "No BAM files specified")));
    }

    // Glob the input bam files
    vector<path> transformedBams = globFiles(bamFiles);

    auto_cpu_timer timer(1, "\nPortcullis prep completed.\nTotal runtime: %ws\n\n");        

    cout << "Running portcullis in prepare mode" << endl
         << "----------------------------------" << endl << endl;

    // Create the prepare class
    Prepare prep(outputDir, strandednessFromString(strandSpecific), force, !copy, useCsi, threads, verbose);
    
    // Prep the input to produce a usable indexed and sorted bam plus, indexed
    // genome and queryable coverage information
    prep.prepare(transformedBams, genomeFile);

    return 0;        
}
    