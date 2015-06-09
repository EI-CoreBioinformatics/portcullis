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

#include <glob.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::vector;

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

#include "samtools_helper.hpp"
#include "genome_mapper.hpp"
#include "portcullis_fs.hpp"
using portcullis::SamtoolsHelper;
using portcullis::PortcullisFS;

#include "prepare.hpp"

bool portcullis::PreparedFiles::valid() {
        
    if (!bfs::exists(getSortedBamFilePath()) && !bfs::symbolic_link_exists(getSortedBamFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find sorted BAM files at: ") + getSortedBamFilePath().string()));
    }

    if (!bfs::exists(getBamIndexFilePath()) && !bfs::symbolic_link_exists(getBamIndexFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find BAM index at: ") + getBamIndexFilePath().string()));
    }

    if (!bfs::exists(getGenomeFilePath()) && !bfs::symbolic_link_exists(getGenomeFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find genome file at: ") + getGenomeFilePath().string()));
    }

    if (!bfs::exists(getGenomeIndexFilePath()) && !bfs::symbolic_link_exists(getGenomeIndexFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find genome index at: ") + getGenomeIndexFilePath().string()));
    }

    if (!bfs::exists(getSettingsFilePath()) && !bfs::symbolic_link_exists(getSettingsFilePath())) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not find settings file at: ") + getSettingsFilePath().string()));
    } 

    return true;
}
    
void portcullis::PreparedFiles::clean() {

    bfs::remove(getUnsortedBamFilePath());
    bfs::remove(getSortedBamFilePath());
    bfs::remove(getBamIndexFilePath());
    bfs::remove(getGenomeFilePath());
    bfs::remove(getGenomeIndexFilePath());
    bfs::remove(getBcfFilePath());
    bfs::remove(getBcfIndexFilePath());
    bfs::remove(getSettingsFilePath());
}
    
portcullis::StrandSpecific portcullis::PreparedFiles::loadSettings() {
    ifstream ifs(getSettingsFilePath().c_str());
    string line;
    while ( std::getline(ifs, line) ) {
        if ( !line.empty() ) {
           if (line.find("SS") == 0) {
               size_t eqPos = line.find("=");                   
               if (eqPos > 0) {
                   string val = (line.substr(eqPos+1)); 
                   boost::trim(val);
                   std::istringstream is(val);
                   string ss;
                   is >> ss;
                   return SSFromString(ss);
               }
           } 
        }
    }
}

    
portcullis::Prepare::Prepare(const path& _outputPrefix, StrandSpecific _strandSpecific, bool _force, bool _useLinks, uint16_t _threads, bool _verbose) {
    output = new PreparedFiles(_outputPrefix);
    strandSpecific = _strandSpecific;
    force = _force;
    useLinks = _useLinks;
    threads = _threads;
    verbose = _verbose;

    if (verbose) {
        cout << "Configured portcullis prep to use the following settings: " << endl
             << " - Output directory: " << output->getPrepDir() << endl
             << " - Strand specific library: " << SSToString(strandSpecific) << endl
             << " - Force prep (cleans output directory): " << boolalpha << force << endl
             << " - Use symbolic links instead of copy where possible: " << boolalpha << useLinks << endl
             << " - Threads (for sorting BAM): " << threads << endl << endl;
    }

    if (force) {
        if (verbose) {
            cout << "Cleaning output dir " << _outputPrefix << " ... ";
            cout.flush();
        }
        output->clean();
        if (verbose) cout << "done." << endl << endl;
    }
}
    
bool portcullis::Prepare::copy(const path& from, const path& to, const string& msg) {

    auto_cpu_timer timer(1, string(" - Copy ") + msg + " - Wall time taken: %ws\n\n");

    bool fileExists = bfs::exists(to) || bfs::symbolic_link_exists(to);

    if (fileExists) {
        if (verbose) cout << "Prepped " << msg << " file detected: " << to << endl;            
    }
    else {

        if (useLinks) {
            bfs::create_symlink(bfs::canonical(from), to);            
            if (verbose) cout << "Created symlink from " << from << " to " << to << endl;
        }
        else {
            if (verbose) {
                cout << "Copying from " << from << " to " << to << " ... ";
                cout.flush();
            }
            std::ifstream  src(from.string(), std::ios::binary);
            std::ofstream  dst(to.string(),   std::ios::binary);

            dst << src.rdbuf();

            src.close();
            dst.close();

            if (verbose) cout << "done." << endl;
        }
    }

    return bfs::exists(to) || bfs::symbolic_link_exists(to);
}
    
bool portcullis::Prepare::genomeIndex() {

    auto_cpu_timer timer(1, " - Genome Index - Wall time taken: %ws\n\n");

    const path genomeFile = output->getGenomeFilePath();
    const path indexFile = output->getGenomeIndexFilePath();

    bool indexExists = bfs::exists(indexFile);

    if (indexExists) {
        if (verbose) cout << "Pre-indexed genome detected: " << indexFile << endl;            
    }
    else {

        if (verbose) {
            cout << "Indexing genome " << genomeFile << " ... ";
            cout.flush();
        }

        // Create the index
        GenomeMapper(genomeFile).buildFastaIndex();

        if (verbose) {
            cout << "done." << endl
                 << "Genome index file created at: " << output->getGenomeIndexFilePath() << endl;                
        }
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

    auto_cpu_timer timer(1, " - BAM Merge - Wall time taken: %ws\n\n");

    path mergedBam = output->getUnsortedBamFilePath();        

    bool mergedBamExists = bfs::exists(mergedBam) || bfs::symbolic_link_exists(mergedBam);

    if (mergedBamExists) {
        if (verbose) cout << "Pre-merged BAM detected: " << mergedBam << endl;
    }
    else {

        if (verbose) {
            cout << "Found " << bamFiles.size() << " BAM files." << endl;
        }

        string mergeCmd = SamtoolsHelper::createMergeBamCmd(bamFiles, mergedBam, threads);

        if (verbose) {
            cout << "Merging BAM using command \"" << mergeCmd << "\" ... ";
            cout.flush();
        }

        int exitCode = system(mergeCmd.c_str());

        if (exitCode != 0 || !bfs::exists(mergedBam)) {
            BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Failed to successfully merge: ") + mergedBam.string()));
        }

        if (verbose) {
        cout << "done." << endl
             << "Merged BAM file created at: " << mergedBam << endl;
        }
    }

    // Return true if the merged BAM exists now, which is should do
    return bfs::exists(mergedBam) || bfs::symbolic_link_exists(mergedBam);        
}

/**
 * Sorts the unsorted bam file if required or forced
 * @param inputBam
 * @return 
 */
bool portcullis::Prepare::bamSort() {

    auto_cpu_timer timer(1, " - BAM Sort - Wall time taken: %ws\n\n");

    const path unsortedBam = output->getUnsortedBamFilePath();
    const path sortedBam = output->getSortedBamFilePath();

    bool sortedBamExists = bfs::exists(sortedBam) || bfs::symbolic_link_exists(sortedBam);

    if (sortedBamExists) {            
        if (verbose) cout << "Pre-sorted BAM detected: " << sortedBam << endl;
    }
    else {

        if (SamtoolsHelper::isCoordSortedBam(unsortedBam.string()) && !force) {

            if (verbose) cout << "Provided BAM appears to be sorted already, just creating symlink instead." << endl;
            bfs::create_symlink(bfs::canonical(unsortedBam), sortedBam);            
            if (verbose) cout << "Created symlink from " << bfs::canonical(unsortedBam) << " to " << sortedBam << endl;
        }
        else {

            // Sort the BAM file by coordinate
            string sortCmd = SamtoolsHelper::createSortBamCmd(unsortedBam, sortedBam, false, threads, "1G");

            if (verbose) {
                cout << "Sorting BAM using command \"" << sortCmd << "\" ... ";
                cout.flush();
            }

            int exitCode = system(sortCmd.c_str());

            path badNameMergeFile = path(sortedBam.string() + ".bam");

            if (bfs::exists(badNameMergeFile) || bfs::symbolic_link_exists(badNameMergeFile)) {
                boost::filesystem::rename(badNameMergeFile, sortedBam);
            }

            if (exitCode != 0 || !bfs::exists(sortedBam) || !SamtoolsHelper::isCoordSortedBam(sortedBam)) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Failed to successfully sort: ") + unsortedBam.string()));
            }

            if (verbose) {
            cout << "done." << endl
                 << "Sorted BAM file created at: " << sortedBam << endl;
            }
        }

    }

    // Return true if the sorted BAM exists now, which is should do
    return bfs::exists(sortedBam) || bfs::symbolic_link_exists(sortedBam);        
}
    
bool portcullis::Prepare::bamIndex() {

    auto_cpu_timer timer(1, " - BAM Index - Wall time taken: %ws\n\n");

    const path sortedBam = output->getSortedBamFilePath();
    const path indexedFile = output->getBamIndexFilePath();

    bool indexedBamExists = bfs::exists(indexedFile) || bfs::symbolic_link_exists(indexedFile);

    if (indexedBamExists) {
        if (verbose) cout << "Pre-indexed BAM detected: " << indexedFile << endl;
    }
    else {

        // Create BAM index
        string indexCmd = SamtoolsHelper::createIndexBamCmd(sortedBam);                

        if (verbose) {
            cout << "Indexing BAM using command \"" << indexCmd << "\" ... ";
            cout.flush();
        }

        int exitCode = system(indexCmd.c_str());                    

        if (exitCode != 0 || !exists(output->getBamIndexFilePath())) {
                BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                        "Failed to successfully index: ") + sortedBam.string()));
        }

        if (verbose) {
            cout << "done." << endl
                 << "BAM index created at: " << output->getBamIndexFilePath() << endl;
        }
    }

    return bfs::exists(indexedFile) || bfs::symbolic_link_exists(indexedFile);
}
    
void portcullis::Prepare::prepare(vector<path> bamFiles, const path& originalGenomeFile) {

    const bool doMerge = bamFiles.size() > 1;        

    if (bamFiles.empty()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "No BAM files to process")));
    }

    path mergedBamFile = doMerge ? output->getUnsortedBamFilePath() : bamFiles[0];

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
                        "Could not copy/symlink BAM file to: ") + output->getUnsortedBamFilePath().string()));
        }            
    }

    // Sort the file
    if (!bamSort()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not sort: ") + output->getUnsortedBamFilePath().string()));
    }

    // Index the sorted file
    if (!bamIndex()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Failed to index: ") + output->getSortedBamFilePath().string()));
    }

    // Copy / Symlink the file to the output dir
    if (!copy(originalGenomeFile, output->getGenomeFilePath(), "genome")) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not copy/symlink genome file to: ") + output->getGenomeFilePath().string()));
    }

    // Test if index exists
    if (!genomeIndex()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                    "Could not create genome map")));
    }

}
   

bool portcullis::Prepare::outputDetails() {

    const path settingsFile = output->getSettingsFilePath();

    std::ofstream outfile(settingsFile.c_str());
    if(!outfile.is_open()) {
        BOOST_THROW_EXCEPTION(PrepareException() << PrepareErrorInfo(string(
                "Could not open settings file for writing: ") + settingsFile.string()));
    }

    outfile << "SS=" << SSToString(strandSpecific) << endl;
    outfile.close();
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
    
int portcullis::Prepare::main(int argc, char *argv[]) {

    // Portcullis args
    vector<path> bamFiles;
    path genomeFile;
    path outputDir;
    bool force;
    string strandSpecific;
    bool useLinks;
    uint16_t threads;
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(helpMessage());
    generic_options.add_options()
            ("output,o", po::value<path>(&outputDir)->default_value(DEFAULT_PREP_OUTPUT_DIR), 
                (string("Output directory for prepared files. Default: ") + DEFAULT_PREP_OUTPUT_DIR).c_str())
            ("force,f", po::bool_switch(&force)->default_value(false), 
                "Whether or not to clean the output directory before processing, thereby forcing full preparation of the genome and bam files.  By default portcullis will only do what it thinks it needs to.")
            ("strand_specific,ss", po::value<string>(&strandSpecific)->default_value(SSToString(StrandSpecific::UNSTRANDED)), 
                "Whether BAM alignments were generated using a strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads)  Default: \"unstranded\"")
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
            ("bam-files,i", po::value< vector<path> >(&bamFiles), "Path to the BAM files to process.")
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
    Prepare prep(outputDir, SSFromString(strandSpecific), force, useLinks, threads, verbose);
    
    // Prep the input to produce a usable indexed and sorted bam plus, indexed
    // genome and queryable coverage information
    prep.prepare(transformedBams, genomeFile);

    // Output any remaining details
    prep.outputDetails();

    return 0;        
}
    