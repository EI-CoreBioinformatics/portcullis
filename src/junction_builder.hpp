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

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
using std::boolalpha;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::mutex;
using std::ofstream;
using std::queue;
using std::thread;
using std::condition_variable;

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

#include <portcullis/intron.hpp>
#include <portcullis/junction.hpp>
#include <portcullis/junction_system.hpp>
using portcullis::Intron;
using portcullis::Junction;
using portcullis::JunctionSystem;

#include "prepare.hpp"
using portcullis::PreparedFiles;



namespace portcullis {

const string DEFAULT_JUNC_OUTPUT = "portcullis_junc/portcullis";
const string DEFAULT_JUNC_SOURCE = "portcullis";
const uint16_t DEFAULT_JUNC_THREADS = 1;

typedef boost::error_info<struct JunctionBuilderError,string> JunctionBuilderErrorInfo;
struct JunctionBuilderException: virtual boost::exception, virtual std::exception { };

struct RegionResult {
    uint64_t splicedCount = 0;
    uint64_t unsplicedCount = 0;
    uint64_t sumQueryLengths = 0;
    int32_t minQueryLength = 100000;
    int32_t maxQueryLength = 0;
    JunctionSystem js;
};

class JunctionBuilder {
private:

    // Can set these from the outside via the constructor
    PreparedFiles prepData;
    path outputDir;
    string outputPrefix;
    uint16_t threads;
    Strandedness strandSpecific;
    bool extra;
    bool separate;
    bool useCsi;
    bool outputExonGFF;
    bool outputIntronGFF;
    string source;
    bool verbose;
    
    // The set of distinct junctions found in the BAM file
    JunctionSystem junctionSystem;
    SplicedAlignmentMap splicedAlignmentMap;
    
    // List of reference sequences (might be shared amongst various objects)
    shared_ptr<RefSeqPtrList> refs;
    
    // Map of reference sequence indicies to reference sequences
    shared_ptr<RefSeqPtrIndexMap> refMap;
    
    // Results from threads
    vector<RegionResult> results;
    
    

protected:
    
    path getUnsplicedBamFile() {
        return path(outputDir.string() + "/" + outputPrefix + ".unspliced.bam");
    }
    
    path getSplicedBamFile() {
        return path(outputDir.string() + "/" + outputPrefix + ".spliced.bam");
    }
    
    path getUnmappedBamFile() {
        return path(outputDir.string() + "/" + outputPrefix + ".unmapped.bam");
    }
    
    path getAssociatedIndexFile(path bamFile) {
        return path(bamFile.string() + ".bai");
    }
        
    void separateBams();
    
    void findJunctions();
    
    void calcExtraMetrics();
        

public:

    
    JunctionBuilder(const path& _prepDir, const path& _output);
    
    virtual ~JunctionBuilder();
    
    string getRefName(const int32_t seqId) { return refs->at(seqId)->name; }
    
    void findJuncs(BamReader& reader, GenomeMapper& gmap, const int32_t seq);
    
    PreparedFiles& getPreparedFiles() { return prepData; }
    
    bool isExtra() const {
        return extra;
    }

    void setExtra(bool extra) {
        this->extra = extra;
    }

    uint16_t getThreads() const {
        return threads;
    }

    void setThreads(uint16_t threads) {
        this->threads = threads;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }
    
    bool isSeparate() const {
        return separate;
    }

    void setSeparate(bool separate) {
        this->separate = separate;
    }

    string getSource() const {
        return source;
    }

    void setSource(string source) {
        this->source = source;
    }
    
    Strandedness getStrandSpecific() const {
        return strandSpecific;
    }

    void setStrandSpecific(Strandedness strandSpecific) {
        this->strandSpecific = strandSpecific;
    }

    bool isUseCsi() const {
        return useCsi;
    }

    void setUseCsi(bool useCsi) {
        this->useCsi = useCsi;
    }
            
    bool isOutputExonGFF() const {
        return outputExonGFF;
    }

    void setOutputExonGFF(bool outputExonGFF) {
        this->outputExonGFF = outputExonGFF;
    }

    bool isOutputIntronGFF() const {
        return outputIntronGFF;
    }

    void setOutputIntronGFF(bool outputIntronGFF) {
        this->outputIntronGFF = outputIntronGFF;
    }


    
    
    /**
     * Populates the set of distinct junctions.  
     * 
     * Also outputs all the unspliced alignments to a separate file if requested
     */
    void process();
    
    static string title() {
        return string("Portcullis Junction Builder Mode Help");
    }
    
    static string description() {
        return string("Analyses all potential junctions found in the input BAM file.\n") +
                      "Run \"portcullis prep ...\" to generate data suitable for junction finding\n" +
                      "before running \"portcullis junc ...\"";
    }
    
    static string usage() {
        return string("portcullis junc [options] <prep_data_dir>");
    }
    
    static int main(int argc, char *argv[]);
};

class JBThreadPool {
public:

    // Constructor.
    JBThreadPool(JunctionBuilder* jb, const uint16_t threads);

    // Destructor.
    ~JBThreadPool();

    // Adds task to a task queue.
    void enqueue(const int32_t index);

    // Shut down the pool.
    void shutDown();

private:
    
    // JunctionBuider
    JunctionBuilder* junctionBuilder;
    
    // Thread pool storage.
    vector<thread> threadPool;

    // Queue to keep track of incoming tasks and task index.
    queue<int32_t> tasks;

    // Task queue mutex.
    mutex tasksMutex;

    // Condition variable.
    condition_variable condition;

    // Indicates that pool needs to be shut down.
    bool terminate;

    // Indicates that pool has been terminated.
    bool stopped;

    // Function that will be invoked by our threads.
    void invoke();
};

}

