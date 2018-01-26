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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <mutex>
using std::boolalpha;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::unique_lock;

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

#include <portcullis/bam/bam_master.hpp>
#include <portcullis/bam/bam_reader.hpp>
#include <portcullis/bam/bam_writer.hpp>
#include <portcullis/bam/depth_parser.hpp>
#include <portcullis/bam/genome_mapper.hpp>
#include <portcullis/seq_utils.hpp>
using namespace portcullis::bam;

#include <portcullis/intron.hpp>
#include <portcullis/junction.hpp>
#include <portcullis/junction_system.hpp>
using portcullis::Intron;
using portcullis::Junction;
using portcullis::JunctionSystem;

#include "junction_builder.hpp"
using portcullis::JBThreadPool;

portcullis::JunctionBuilder::JunctionBuilder(const path& _prepDir, const path& _output) {
	prepData = PreparedFiles(_prepDir);
	outputDir = _output.empty() ? path(".") : _output.parent_path();
	outputPrefix = _output.empty() ? "portcullis" : _output.leaf().string();
	threads = 1;
	extra = false;
	useCsi = false;
	strandSpecific = Strandedness::UNKNOWN;
	source = "portcullis";
	verbose = false;
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
	// Ensure output directory exists
	if (!bfs::exists(outputDir)) {
		if (!bfs::create_directories(outputDir)) {
			BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
									  "Could not create output directory at: ") + outputDir.string()));
		}
	}
	if (!bfs::exists(prepData.getSortedBamFilePath())) {
		BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
								  "Could not find prepared BAM file at: ") + prepData.getSortedBamFilePath().string()));
	}
	// Test if we have all the required data
	if (!prepData.valid(useCsi)) {
		BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
								  "Prepared data is not complete: ") + prepData.getPrepDir().string()));
	}
	const path sortedBamFile = prepData.getSortedBamFilePath();
	// Acquire list of reference sequences
	BamReader reader(prepData.getSortedBamFilePath());
	reader.open();
	refs = reader.createRefList();
	refMap = reader.createRefMap(*refs);
	reader.close();
	junctionSystem.setRefs(refs);
	if (refs->size() < threads) {
		cerr << "Warning: User requested " << threads << " threads but there are only " << refs->size() << " target sequences to process.  Setting number of threads to " << refs->size() << "." << endl << endl;
		threads = refs->size();
	}
	// Must separate BAMs if extra metrics are requested
	if (extra && !separate) {
		separate = true;
		cerr << "Warning: User requested that separated BAMS should not be output but user did request extra metrics to be calculated.  This requires separated BAMs to be produced." << endl << endl;
	}
	// Output settings requested
	cout << "Settings:" << endl
		 << std::boolalpha
		 << " - BAM Strandedness: " << strandednessToString(strandSpecific) << endl
		 << " - BAM Read Orientation: " << orientationToString(orientation) << endl
		 << " - BAM Indexing mode: " << (useCsi ? "CSI" : "BAI") << endl
		 << " - Threads: " << threads << endl
		 << " - Separate BAMs: " << separate << endl
		 //<< " - Calculate additional metrics: " << extra << endl
		 << endl;
	cout << reader.bamDetails() << endl;
	// Separate spliced from unspliced reads and save to file if requested
	if (separate) {
		separateBams();
	}
	// The core interesting work is done here
	findJunctions();
	if (extra) {
		calcExtraMetrics();
	}
	cout << "Saving junctions: " << endl;
	junctionSystem.saveAll(path(outputDir.string() + "/" + outputPrefix), source, false, this->outputExonGFF, this->outputIntronGFF);

	// Also do a strand analysis as this is cheap and quick to do.
	std::pair<Orientation, Strandedness> actual_config = junctionSystem.determineStrandedness(true);
	Orientation actual_orientation = actual_config.first;
	Strandedness actual_strandedness = actual_config.second;
	cout << "Determined sequence orientation to be: " << orientationToLongString(actual_orientation) << endl;
	cout << "Determined RNAseq strandedness to be: " << strandednessToLongString(actual_strandedness) << endl << endl;
	if (strandSpecific != Strandedness::UNKNOWN && strandSpecific != actual_strandedness) {
		cerr << "Warning!  User input and portcullis disagree about the strandedness of the dataset" << endl << endl;
	}
}

void portcullis::JunctionBuilder::separateBams() {
	auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
	uint64_t splicedCount = 0;
	uint64_t unsplicedCount = 0;
	uint64_t unmappedCount = 0;
	const path unsplicedFile = getUnsplicedBamFile();
	const path splicedFile = getSplicedBamFile();
	const path unmappedFile = getUnmappedBamFile();
	BamWriter unsplicedWriter(unsplicedFile);
	BamWriter splicedWriter(splicedFile);
	BamWriter unmappedWriter(unmappedFile);
	BamReader reader(prepData.getSortedBamFilePath());
	reader.open();
	cout << "Splitting BAM:" << endl;
	if (threads > 1) {
		cerr << " - Note: this part of the pipeline is single threaded." << endl;
	}
	cout << " - Saving unspliced alignments to: " << unsplicedFile << endl;
	unsplicedWriter.open(reader.getHeader());
	cout << " - Saving spliced alignments to: " << splicedFile << endl;
	splicedWriter.open(reader.getHeader());
	cout << " - Saving unmapped reads to: " << unmappedFile << endl;
	unmappedWriter.open(reader.getHeader());
	cout << " - Processing BAM ...";
	cout.flush();
	while (reader.next()) {
		const BamAlignment& al = reader.current();
		if (al.isSplicedRead()) {
			splicedWriter.write(al);
			splicedCount++;
			if (extra) {
				// Record alignment name in map
				size_t code = std::hash<string>()(al.deriveName());
				splicedAlignmentMap[code]++;
			}
		}
		else if (al.isMapped()) {
			unsplicedWriter.write(al);
			unsplicedCount++;
		}
		else {
			unmappedWriter.write(al);
			unmappedCount++;
		}
	}
	cout << " done." << endl;
	cout << " - Found " << splicedCount << " spliced alignments." << endl;
	cout << " - Found " << unsplicedCount << " unspliced alignments." << endl;
	cout << " - Found " << unmappedCount << " unmapped reads." << endl;
	reader.close();
	unsplicedWriter.close();
	splicedWriter.close();
	unmappedWriter.close();
	cout << " - Indexing unspliced alignments ... ";
	cout.flush();
	// Create BAM index
	string unsplicedIndexCmd = BamHelper::createIndexBamCmd(unsplicedFile, useCsi);
	int unsExitCode = system(unsplicedIndexCmd.c_str());
	if (unsExitCode != 0) {
        BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
									  "Problem indexing unspliced BAM: ") + unsplicedFile.string()));
    }

    cout << "done." << endl;
	cout << " - Indexing spliced alignments ... ";
	cout.flush();
	// Create BAM index
	string splicedIndexCmd = BamHelper::createIndexBamCmd(splicedFile, useCsi);
	int sExitCode = system(splicedIndexCmd.c_str());
    if (sExitCode != 0) {
        BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
									  "Problem indexing spliced BAM: ") + splicedFile.string()));
    }
	cout << "done." << endl;
}

void portcullis::JunctionBuilder::findJunctions() {
	auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
	// Add each target sequence as a chunk of work for the thread pool
	results.clear();
	results.resize(refs->size());
	// Create the thread pool and start the threads
	cout << "Creating " << threads << " threads, each with BAM and genome indicies loaded ...";
	cout.flush();
	JBThreadPool pool(this, threads);
	cout << " done." << endl;
	cout << "Finding junctions and calculating basic metrics:" << endl;
	cout << " - Queueing " << refs->size() << " target sequences for processing in the thread pool" << endl;
	cout << " - Processing: " << endl;
	for (size_t i = 0; i < refs->size(); i++) {
		results[i].js.setRefs(refs); // Make sure junction system has reference sequence list available
		results[i].name = refs->at(i)->name;
		pool.enqueue(refs->at(i)->index);
	}
	// Waits for all threads to complete
	pool.shutDown();
	cout << " - All threads completed." << endl << " - Combining results from threads." << endl << endl;
	uint64_t unsplicedCount = 0;
	uint64_t splicedCount = 0;
	uint64_t sumQueryLengths = 0;
	int32_t minQueryLength = INT32_MAX;
	int32_t maxQueryLength = 0;
	cout << std::left << std::setw(12) << "Sequence" << "\t"
		 << std::right << std::setw(12) << "unspliced" << "\t"
		 << std::right << std::setw(12) << "spliced" << "\t"
		 << std::right << std::setw(12) << "total" << endl;
	for (auto & res : results) {
		junctionSystem.append(res.js);
		unsplicedCount += res.unsplicedCount;
		splicedCount += res.splicedCount;
		sumQueryLengths += res.sumQueryLengths;
		minQueryLength = min(minQueryLength, res.minQueryLength);
		maxQueryLength = max(maxQueryLength, res.maxQueryLength);
		cout << std::left << std::setw(12) << res.name << "\t"
			 << std::right << std::setw(12) << res.unsplicedCount << "\t"
			 << std::right << std::setw(12) << res.splicedCount << "\t"
			 << std::right << std::setw(12) << res.splicedCount + res.unsplicedCount << endl;
	}
	cout << endl << "Sorting and reindexing merged junctions...";
	cout.flush();
	junctionSystem.sort(); // Make sure the output is properly ordered
	junctionSystem.index(); // Add unique identifiers to each junction
	cout << " done." << endl << endl;
	// Calculate some alignment stats
	uint64_t totalAlignments = splicedCount + unsplicedCount;
	double meanQueryLength = (double) sumQueryLengths / (double) totalAlignments;
	junctionSystem.setQueryLengthStats(minQueryLength, meanQueryLength, maxQueryLength);
	cout << "Final stats:" << endl
		 << " - Processed " << totalAlignments << " alignments." << endl
		 << " - Alignment query length statistics: min: " << minQueryLength << "; mean: " << meanQueryLength << "; max: " << maxQueryLength << ";" << endl
		 << " - Found " << junctionSystem.size() << " junctions from " << splicedCount << " spliced alignments." << endl
		 << " - Found " << unsplicedCount << " unspliced alignments." << endl;
	// Calculate additional junction stats
	if (junctionSystem.size() > 1) {
		cout << " - Calculating junctions stats that require comparisons with other junctions...";
		cout.flush();
		junctionSystem.calcJunctionStats();
		cout << " done." << endl;
	}
}

void portcullis::JunctionBuilder::calcExtraMetrics() {
	auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
	cout << "Calculating extra junction metrics:" << endl;
	if (threads > 1) {
		cerr << " - Note: this part of the pipeline is single threaded." << endl;
	}
	// Requires BAMs to be separated
	cout << " - Calculating multiple mapping stats ...";
	cout.flush();
	junctionSystem.calcMultipleMappingStats(splicedAlignmentMap);
	cout << " done" << endl;
	// Count the number of alignments found in upstream and downstream flanking
	// regions for each junction
	cout << " - Analysing unspliced alignments around junctions ...";
	cout.flush();
	junctionSystem.findFlankingAlignments(getUnsplicedBamFile());
	cout << " - Calculating unspliced alignment coverage around junctions ...";
	cout.flush();
	junctionSystem.calcCoverage(getUnsplicedBamFile(), strandSpecific);
}

void portcullis::JunctionBuilder::findJuncs(BamReader& reader, GenomeMapper& gmap, int32_t seq) {
	uint64_t splicedCount = 0;
	uint64_t unsplicedCount = 0;
	uint32_t lastCalculatedJunctionIndex = 0;
	uint64_t sumQueryLengths = 0;
	int32_t minQueryLength = INT32_MAX;
	int32_t maxQueryLength = 0;
	reader.setRegion(seq, 0, refs->at(seq)->length);
	while (reader.next()) {
		const BamAlignment& al = reader.current();
		while (results[seq].js.size() > 0 && lastCalculatedJunctionIndex < results[seq].js.size() &&
				al.getPosition() > results[seq].js.getJunctionAt(lastCalculatedJunctionIndex)->getIntron()->end) {
			JunctionPtr j = results[seq].js.getJunctionAt(lastCalculatedJunctionIndex);
			j->calcMetrics(this->orientation);
			j->processJunctionWindow(gmap);
			j->clearAlignments();
			lastCalculatedJunctionIndex++;
		}
		// Calc alignment stats
		int32_t len = al.getLength();
		minQueryLength = min(minQueryLength, len);
		maxQueryLength = max(maxQueryLength, len);
		sumQueryLengths += len;
		if (results[seq].js.addJunctions(al)) {
			splicedCount++;
		}
		else {
			unsplicedCount++;
		}
	}
	while (results[seq].js.size() > 0 && lastCalculatedJunctionIndex < results[seq].js.size()) {
		JunctionPtr j = results[seq].js.getJunctionAt(lastCalculatedJunctionIndex);
		j->calcMetrics(this->orientation);
		j->processJunctionWindow(gmap);
		j->clearAlignments();
		lastCalculatedJunctionIndex++;
	}
	// Update result vector
	results[seq].splicedCount = splicedCount;
	results[seq].unsplicedCount = unsplicedCount;
	results[seq].minQueryLength = minQueryLength;
	results[seq].maxQueryLength = maxQueryLength;
	results[seq].sumQueryLengths = sumQueryLengths;
}

int portcullis::JunctionBuilder::main(int argc, char *argv[]) {
	// Portcullis args
	string prepDir;
	string output;
	uint16_t threads;
	bool extra;
	bool separate;
	string strandSpecific;
	string orientation;
	bool useCsi;
	bool exongff;
	bool introngff;
	string source;
	bool verbose;
	bool help;
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	// Declare the supported options.
	po::options_description system_options("System options", w.ws_col, (unsigned) ((double) w.ws_col / 1.5));
	system_options.add_options()
	("threads,t", po::value<uint16_t>(&threads)->default_value(1),
	 "The number of threads to use.  Note that increasing the number of threads will also increase memory requirements.")
	("separate", po::bool_switch(&separate)->default_value(false),
	 "Separate spliced from unspliced reads.  Creates two new BAM files.")
	("orientation", po::value<string>(&orientation)->default_value(orientationToString(Orientation::UNKNOWN)),
	 "The orientation of the reads that produced the BAM alignments: \"F\" (Single-end forward orientation); \"R\" (single-end reverse orientation); \"FR\" (paired-end, with reads sequenced towards center of fragment -> <-.  This is usual setting for most Illumina paired end sequencing); \"RF\" (paired-end, reads sequenced away from center of fragment <- ->); \"FF\" (paired-end, reads both sequenced in forward orientation); \"RR\" (paired-end, reads both sequenced in reverse orientation); \"UNKNOWN\" (default, portcullis will workaround any calculations requiring orientation information)")
	("strandedness", po::value<string>(&strandSpecific)->default_value(strandednessToString(Strandedness::UNKNOWN)),
	 "Whether BAM alignments were generated using a type of strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads); \"UNKNOWN\" (default, portcullis will workaround any calculations requiring strandedness information)")
	("use_csi,c", po::bool_switch(&useCsi)->default_value(false),
	 "Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it supports very long target sequences (probably not an issue unless you are working on huge genomes).  BAI has the advantage that it is more widely supported (useful for viewing in genome browsers).")
	("verbose,v", po::bool_switch(&verbose)->default_value(false),
	 "Print extra information")
	("help", po::bool_switch(&help)->default_value(false), "Produce help message")
	;
	po::options_description output_options("Output options", w.ws_col, (unsigned) ((double) w.ws_col / 1.5));
	output_options.add_options()
	("output,o", po::value<string>(&output)->default_value(DEFAULT_JUNC_OUTPUT),
	 "Output prefix for files generated by this program.")
	("exon_gff", po::bool_switch(&exongff)->default_value(false),
	 "Output exon-based junctions in GFF format.")
	("intron_gff", po::bool_switch(&introngff)->default_value(false),
	 "Output intron-based junctions in GFF format.")
	("source", po::value<string>(&source)->default_value(DEFAULT_JUNC_SOURCE),
	 "The value to enter into the \"source\" field in GFF files.")
	;
	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	po::options_description hidden_options("Hidden options");
	hidden_options.add_options()
	("extra", po::bool_switch(&extra)->default_value(false),
	 "Calculate additional metrics that take some time to generate.  Automatically activates BAM splitting mode (--separate).")
	("prep_data_dir,i", po::value<string>(&prepDir), "Path to directory containing prepared data.")
	;
	// Positional option for the input bam file
	po::positional_options_description p;
	p.add("prep_data_dir", 1);
	// Combine non-positional options
	po::options_description display_options;
	display_options.add(system_options).add(output_options);
	// Combine non-positional options
	po::options_description cmdline_options;
	cmdline_options.add(display_options).add(hidden_options);
	// Parse command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);
	// Output help information the exit if requested
	if (help || argc <= 1) {
		cout << title() << endl << endl
			 << description() << endl << endl
			 << "Usage: " << usage() << endl
			 << display_options << endl;
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
	JunctionBuilder jb(prepDir, output);
	jb.setThreads(threads);
	jb.setExtra(extra);
	jb.setSeparate(separate);
	jb.setSource(source);
	jb.setStrandSpecific(strandednessFromString(strandSpecific));
	jb.setOrientation(orientationFromString(orientation));
	jb.setUseCsi(useCsi);
	jb.setOutputExonGFF(exongff);
	jb.setOutputIntronGFF(introngff);
	jb.setVerbose(verbose);
	jb.process();
	return 0;
}


// ********* Thread Pool ************

portcullis::JBThreadPool::JBThreadPool(JunctionBuilder* jb, const uint16_t threads) : terminate(false), stopped(false) {
	junctionBuilder = jb;
	// Create number of required threads and add them to the thread pool vector.
	for (int i = 0; i < threads; i++) {
		// Add the thread onto the thread pool (providing an index so we can get the the correct BAM reader again)
		threadPool.emplace_back(thread(&portcullis::JBThreadPool::invoke, this));
	}
}

void portcullis::JBThreadPool::enqueue(const int32_t index) {
	// Scope based locking.
	{
		// Put unique lock on task mutex.
		unique_lock<mutex> lock(tasksMutex);
		// Push task into queue.
		tasks.push(index);
	}
	// Wake up one thread.
	condition.notify_one();
}

void portcullis::JBThreadPool::invoke() {
	// Create the genome mapper
	GenomeMapper gmap(junctionBuilder->getPreparedFiles().getGenomeFilePath());
	// Load the fasta index
	gmap.loadFastaIndex();
	// Create a BAM reader for this thread
	BamReader reader(junctionBuilder->getPreparedFiles().getSortedBamFilePath());
	// Open the BAM file... this will load the index, which might take some time on large BAMs
	reader.open();
	int32_t id;
	while (true) {
		// Scope based locking.
		{
			// Put unique lock on task mutex.
			unique_lock<mutex> lock(tasksMutex);
			// Wait until queue is not empty or termination signal is sent.
			condition.wait(lock, [this] {
				return !tasks.empty() || terminate;
			});
			// If termination signal received and queue is empty then exit else continue clearing the queue.
			// Make sure we close the reader before exiting
			if (terminate && tasks.empty()) {
				reader.close();
				return;
			}
			// Get next task in the queue.
			id = tasks.front();
			cout << "   - " << junctionBuilder->getRefName(id) << endl;
			// Remove it from the queue.
			tasks.pop();
		}
		// Execute the task.
		junctionBuilder->findJuncs(reader, gmap, id);
	}
}

void portcullis::JBThreadPool::shutDown() {
	// Scope based locking.
	{
		// Put unique lock on task mutex.
		unique_lock<mutex> lock(tasksMutex);
		// Set termination flag to true.
		terminate = true;
	}
	// Wake up all threads.
	condition.notify_all();
	// Join all threads.
	for (thread & thread : threadPool) {
		thread.join();
	}
	// Empty workers vector.
	threadPool.empty();
	// Indicate that the pool has been shut down.
	stopped = true;
}

// Destructor.

portcullis::JBThreadPool::~JBThreadPool() {
	if (!stopped) {
		shutDown();
	}
}
