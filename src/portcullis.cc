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

#include <execinfo.h>
#include <signal.h>
#include <sys/ioctl.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
using boost::to_upper_copy;
using boost::timer::auto_cpu_timer;
namespace po = boost::program_options;

#include <portcullis/portcullis_fs.hpp>
using portcullis::PortcullisFS;

#include "junction_builder.hpp"
#include "prepare.hpp"
#include "junction_filter.hpp"
#include "bam_filter.hpp"
using portcullis::JunctionBuilder;
using portcullis::Prepare;
using portcullis::JunctionFilter;
using portcullis::BamFilter;

typedef boost::error_info<struct PortcullisError, string> PortcullisErrorInfo;
struct PortcullisException: virtual boost::exception, virtual std::exception { };

// Default values for arguments
const uint16_t DEFAULT_THREADS = 4;
const uint32_t DEFAULT_CHUNK_SIZE_PER_THREAD = 10000;
const uint32_t DEFAULT_GAP_SIZE = 100;

// Global variable! :(
portcullis::PortcullisFS portcullis::pfs;

enum class Mode {
	PREP,
	JUNC,
	FILTER,
	BAM_FILT,
    FULL
};

void print_backtrace(int depth = 0) {
    try {
        throw;
    }    
    // this block shows how to handle exceptions of some known type
    // You can have your own types instead of std::exception
    catch (const std::exception & e) {
        std::cerr << std::string(depth, ' ') << e.what() << std::endl;
        try {
            std::rethrow_if_nested(e);
        } catch (...) {
            print_backtrace(++depth);
        }
    }    
    // Not all nesting exceptions will be of a known type, but if they use the
    // mixin type std::nested_exception, then we can at least handle them enough to
    // get the nested exception:
    catch (const std::nested_exception & ne) {
        std::cerr << std::string(depth, ' ') << "Unknown nesting exception\n";

        try {
            ne.rethrow_nested();
        } catch (...) {
            print_backtrace(++depth);
        }
    }    
    // Exception nesting works through inheritance, which means that if you
    // can't inherit from the type, then you can't 'mixin' std::nesting exception.
    // If you try something like std::throw_with_nested( int{10} ); Then you'll
    // hit this catch block when printing the backtrace.
    catch (...) {
        std::cerr << std::string(depth, ' ') << "Unknown exception\n";
    }
}

    

Mode parseMode(string mode) {
	string upperMode = boost::to_upper_copy(mode);
	if (upperMode == string("PREP") || upperMode == string("PREPARE")) {
		return Mode::PREP;
	}
	else if (upperMode == string("JUNC") || upperMode == string("ANALYSE") || upperMode == string("ANALYZE")) {
		return Mode::JUNC;
	}
	else if (upperMode == string("FILTER") || upperMode == string("FILT")) {
		return Mode::FILTER;
	}
	else if (upperMode == string("BAMFILT")) {
		return Mode::BAM_FILT;
	}
	else if (upperMode == string("FULL")) {
		return Mode::FULL;
	}
    else {
		BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
								  "Could not recognise mode string: ") + mode));
	}
}

string title() {
	return string("Portcullis Help");
}

string description() {
	return string("Portcullis is a tool to identify genuine splice junctions using aligned RNAseq reads");
}

string usage() {
	return string("portcullis [options] <mode> <mode_args>");
}

string modes() {
	return string(
			   " - full    - Full pipeline.  Runs prep, junc, filt (and optionally bamfilt) in sequence\n") +
		   " - prep    - Step 1: Prepares a genome and bam file(s) ready for junction analysis\n" +
		   " - junc    - Step 2: Perform junction analysis on prepared data\n" +
		   " - filt    - Step 3: Discard unlikely junctions\n" +
		   " - bamfilt - Step 4: Filters a BAM to remove any reads associated with invalid\n" +
		   "             junctions";
}

string fulltitle() {
	return string("Portcullis Full Pipeline Mode Help");
}

string fulldescription() {
	return string("Runs prep, junc, filter, and optionally bamfilt, as a complete pipeline.  Assumes\n") +
		   "that the self-trained machine learning approach for filtering is to be used.";
}

string fullusage() {
	return string("portcullis full [options] <genome-file> (<bam-file>)+");
}


int mainFull(int argc, char *argv[]) {
	// Portcullis args
	std::vector<path> bamFiles;
	path genomeFile;
	path outputDir;
	path referenceFile;
	string strandSpecific;
	string orientation;
	uint16_t threads;
	bool separate;
    bool extra;
    bool copy;
	bool keep_temp;
    bool force;
	bool useCsi;
	bool saveBad;
	bool exongff;
	bool introngff;
	bool bamFilter;
    string source;
	uint32_t max_length;
	uint32_t mincov;
	string canonical;
	bool verbose;
	bool help;
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	// Declare the supported options.
	po::options_description system_options("System options", w.ws_col, (unsigned)((double)w.ws_col / 1.5));
	system_options.add_options()
	("threads,t", po::value<uint16_t>(&threads)->default_value(1),
	 "The number of threads to use.  Note that increasing the number of threads will also increase memory requirements.  Default: 1")
	("verbose,v", po::bool_switch(&verbose)->default_value(false),
	 "Print extra information")
	("help", po::bool_switch(&help)->default_value(false), "Produce help message")
	;
	po::options_description output_options("Output options", w.ws_col, (unsigned)((double)w.ws_col / 1.5));
	output_options.add_options()
	("output,o", po::value<path>(&outputDir)->default_value("portcullis_out"),
	 "Output directory. Default: portcullis_out")
	("bam_filter,b", po::bool_switch(&bamFilter)->default_value(false),
	 "Filter out alignments corresponding with false junctions.  Warning: this is time consuming; make sure you really want to do this first!")
	("exon_gff", po::bool_switch(&exongff)->default_value(false),
	 "Output exon-based junctions in GFF format.")
	("intron_gff", po::bool_switch(&introngff)->default_value(false),
	 "Output intron-based junctions in GFF format.")
	("source", po::value<string>(&source)->default_value("portcullis"),
	 "The value to enter into the \"source\" field in GFF files.")
	;
	po::options_description prepare_options("Input options", w.ws_col, (unsigned)((double)w.ws_col / 1.5));
	prepare_options.add_options()
	("force", po::bool_switch(&force)->default_value(false),
	 "Whether or not to clean the output directory before processing, thereby forcing full preparation of the genome and bam files.  By default portcullis will only do what it thinks it needs to.")
	("copy", po::bool_switch(&copy)->default_value(false),
	 "Whether to copy files from input data to prepared data where possible, otherwise will use symlinks.  Will require more time and disk space to prepare input but is potentially more robust.")
	("keep_temp", po::bool_switch(&keep_temp)->default_value(false),
	 "Whether keep any temporary files created during the prepare stage of portcullis.  This might include BAM files and indexes.")
	("use_csi", po::bool_switch(&useCsi)->default_value(false),
	 "Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it supports very long target sequences (probably not an issue unless you are working on huge genomes).  BAI has the advantage that it is more widely supported (useful for viewing in genome browsers).")
	;
    
    po::options_description analysis_options("Analysis options", w.ws_col, (unsigned)((double)w.ws_col / 1.5));
	analysis_options.add_options()
	("orientation", po::value<string>(&orientation)->default_value(orientationToString(Orientation::UNKNOWN)),
	 "The orientation of the reads that produced the BAM alignments: \"F\" (Single-end forward orientation); \"R\" (single-end reverse orientation); \"FR\" (paired-end, with reads sequenced towards center of fragment -> <-.  This is usual setting for most Illumina paired end sequencing); \"RF\" (paired-end, reads sequenced away from center of fragment <- ->); \"FF\" (paired-end, reads both sequenced in forward orientation); \"RR\" (paired-end, reads both sequenced in reverse orientation); \"UNKNOWN\" (default, portcullis will workaround any calculations requiring orientation information)")
	("strandedness", po::value<string>(&strandSpecific)->default_value(strandednessToString(Strandedness::UNKNOWN)),
	 "Whether BAM alignments were generated using a type of strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads); \"UNKNOWN\" (default, portcullis will workaround any calculations requiring strandedness information)")
	("separate", po::bool_switch(&separate)->default_value(false),
	 "Separate spliced from unspliced reads.")
	("extra", po::bool_switch(&extra)->default_value(false),
	 "Calculate additional metrics that take some time to generate.  Automatically activates BAM splitting mode (--separate).")
    ;
	
	po::options_description filter_options("Filtering options", w.ws_col, (unsigned)((double)w.ws_col / 1.5));
	filter_options.add_options()
	("reference,r", po::value<path>(&referenceFile),
	 "Reference annotation of junctions in BED format.  Any junctions found by the junction analysis tool will be preserved if found in this reference file regardless of any other filtering criteria.  If you need to convert a reference annotation from GTF or GFF to BED format portcullis contains scripts for this.")
	("max_length", po::value<uint32_t>(&max_length)->default_value(0),
	 "Filter junctions longer than this value.  Default (0) is to not filter based on length.")
	("canonical", po::value<string>(&canonical)->default_value("OFF"),
	 "Keep junctions based on their splice site status.  Valid options: OFF,C,S,N. Where C = Canonical junctions (GT-AG), S = Semi-canonical junctions (AT-AC, or  GC-AG), N = Non-canonical.  OFF means, keep all junctions (i.e. don't filter by canonical status).  User can separate options by a comma to keep two categories.")
	("min_cov", po::value<uint32_t>(&mincov)->default_value(1),
	 "Only keep junctions with a number of split reads greater than or equal to this number")
	("save_bad", po::bool_switch(&saveBad)->default_value(false),
	 "Saves bad junctions (i.e. junctions that fail the filter), as well as good junctions (those that pass)")
	;
    
	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	po::options_description hidden_options("Hidden options");
	hidden_options.add_options()
	("bam-files", po::value< std::vector<path> >(&bamFiles), "Path to the BAM files to process.")
	("genome-file", po::value<path>(&genomeFile), "Path to the genome file to process.")
	;
	// Positional option for the input bam file
	po::positional_options_description p;
	p.add("genome-file", 1);
	p.add("bam-files", -1);
	// Combine non-positional options for displaying to the user
	po::options_description display_options;
	display_options.add(system_options).add(output_options).add(prepare_options).add(analysis_options).add(filter_options);
	// Combine non-positional options for use at the command line
	po::options_description cmdline_options;
	cmdline_options.add(display_options).add(hidden_options);
	// Parse command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);
	// Output help information the exit if requested
	if (help || argc <= 1) {
		cout << fulltitle() << endl << endl
			 << fulldescription() << endl << endl
			 << "Usage: " << fullusage() << endl
			 << display_options << endl << endl;
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
		if (!create_directories(outputDir)) {
			BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
									  "Could not create output directory: ") + outputDir.string()));
		}
	}
    
	// ************ Prepare input data (BAMs + genome) ***********
	cout << "Preparing input data (BAMs + genome)" << endl
		 << "----------------------------------" << endl << endl;
	path prepDir = path(outputDir.string() + "/1-prep");
	// Create the prepare class
	Prepare prep(prepDir);
	prep.setForce(force);
	prep.setUseLinks(!copy);
	prep.setUseCsi(useCsi);
	prep.setThreads(threads);
	prep.setVerbose(verbose);
	// Prep the input to produce a usable indexed and sorted bam plus, indexed
	// genome and queryable coverage information
	prep.prepare(transformedBams, genomeFile);
    
	// ************ Identify all junctions and calculate metrics ***********
	cout << "Identifying junctions and calculating metrics" << endl
		 << "---------------------------------------------" << endl << endl;
	path juncDir = outputDir.string() + "/2-junc";
	path juncOut = juncDir.string() + "/portcullis_all";
	// Identify junctions and calculate metrics
	JunctionBuilder jb(prepDir, juncOut);
	jb.setThreads(threads);
	jb.setExtra(false);     // Run in fast mode
	jb.setSeparate(false);  // Run in fast mode
	jb.setStrandSpecific(strandednessFromString(strandSpecific));
	jb.setOrientation(orientationFromString(orientation));
    jb.setExtra(extra);
    jb.setSeparate(separate);
	jb.setSource(source);
	jb.setUseCsi(useCsi);
	jb.setOutputExonGFF(exongff);
	jb.setOutputIntronGFF(introngff);
	jb.setVerbose(verbose);
	jb.process();
    
	// ************ Use default filtering strategy *************
	cout << "Filtering junctions" << endl
		 << "-------------------" << endl << endl;
	path filtOut = outputDir.string() + "/3-filt/portcullis_filtered";
	path juncTab = juncDir.string() + "/portcullis_all.junctions.tab";
	JunctionFilter filter(prepDir, juncTab, filtOut);
	filter.setVerbose(verbose);
	filter.setSource(source);
	filter.setMaxLength(max_length);
	filter.setCanonical(canonical);
	filter.setMinCov(mincov);
	filter.setTrain(true);
	filter.setThreads(threads);
	filter.setENN(false);
	filter.setOutputExonGFF(exongff);
	filter.setOutputIntronGFF(introngff);
	filter.setSaveBad(saveBad);
	filter.filter();
    
	// *********** BAM filter *********
	if (bamFilter) {
		cout << "Filtering BAMs" << endl
			 << "--------------" << endl << endl;
		path filtJuncTab = path(filtOut.string() + ".pass.junctions.tab");
		path bamFile = path(prepDir.string() + "/portcullis.sorted.alignments.bam");
		path filteredBam = path(outputDir.string() + "/portcullis.filtered.bam");
		BamFilter bamFilter(filtJuncTab.string(), bamFile.string(), filteredBam.string());
		//bamFilter.setStrandSpecific(strandednessFromString(strandSpecific));
		//bamFilter.setOrientation(orientationFromString(orientation));
		bamFilter.setUseCsi(useCsi);
		bamFilter.setVerbose(verbose);
		bamFilter.filter();
	}
    
    if (!keep_temp) {
        cout << "Cleaning temporary files...";
        cout.flush();			 
        prep.getOutput()->clean();
        cout << " done." << endl << endl;
    }
    
	return 0;
}


void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    fprintf(stderr, "Stack trace:\n");
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}


/**
 * Start point for portcullis.
 */
int main(int argc, char *argv[]) {
	
    signal(SIGSEGV, handler);
    try {
		// Portcullis args
		string modeStr;
		std::vector<string> others;
		bool verbose;
		bool version;
		bool help;
		struct winsize w;
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
		// Declare the supported options.
		po::options_description generic_options("Options", w.ws_col, (unsigned)((double)w.ws_col / 1.7));
		generic_options.add_options()
		("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print extra information")
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
		// Always output version information but exit if version info was all user requested
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "Portcullis"
#endif
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.X"
#endif
		portcullis::pfs = PortcullisFS(argv[0]);
		portcullis::pfs.setVersion(PACKAGE_VERSION);
		// End if verbose was requested at this level, outputting file system details.
		if (verbose) {
			cout << endl
				 << "Project filesystem" << endl
				 << "------------------" << endl
				 << portcullis::pfs << endl;
			//return 0;
		}
		// Output help information the exit if requested
		if (argc == 1 || (argc == 2 && verbose) || (argc == 2 && help) || (argc == 3 && verbose && help)) {
			cout << title() << endl << endl
				 << description() << endl << endl
				 << "Modes:" << endl << modes() << endl << endl
				 << "Usage: " << usage() << endl << endl
				 << generic_options << endl << endl;
			return 1;
		}
		// End if version was requested.
		if (version) {
			cout << PACKAGE_NAME << " " << PACKAGE_VERSION << endl;
			return 0;
		}
		else {
			cout << "Portcullis V" << PACKAGE_VERSION << endl << endl;
		}
		// If we've got this far parse the command line properly
		Mode mode = parseMode(modeStr);
		const int modeArgC = argc - 1;
		char** modeArgV = argv + 1;
		// Set static variables in downstream subtools so they know where to get their resources from
		JunctionFilter::dataDir = portcullis::pfs.getDataDir();
		JunctionSystem::version = portcullis::pfs.getVersion();
		if (mode == Mode::PREP) {
			Prepare::main(modeArgC, modeArgV);
		}
		else if (mode == Mode::JUNC) {
			JunctionBuilder::main(modeArgC, modeArgV);
		}
		else if (mode == Mode::FILTER) {
			JunctionFilter::main(modeArgC, modeArgV);
		}
		else if (mode == Mode::BAM_FILT) {
			BamFilter::main(modeArgC, modeArgV);
		}
        else if (mode == Mode::FULL) {
			mainFull(modeArgC, modeArgV);
		}
		else {
			BOOST_THROW_EXCEPTION(PortcullisException() << PortcullisErrorInfo(string(
									  "Unrecognised portcullis mode: ") + modeStr));
		}
	}
	catch (po::error& e) {
		cerr << "Error: Parsing Command Line: " << e.what() << endl;
		return 1;
	}
	catch (boost::exception& e) {
		cerr << boost::diagnostic_information(e);
                return 4;
	}
	catch (std::exception& e) {
		cerr << "Error: " << e.what() << endl;
                print_backtrace();
		return 5;
	}
	catch (const char* msg) {
		cerr << "Error: " << msg << endl;
                print_backtrace();
		return 6;
	}
	catch (...) {
		cerr << "Error: Exception of unknown type!" << endl;
                print_backtrace();
		return 7;
	}
	return 0;
}

