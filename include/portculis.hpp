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
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include <boost/unordered_map.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;

#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;

#include "genome_mapper.hpp"
#include "intron.hpp"
#include "junction.hpp"
#include "junction_system.hpp"
using portculis::GenomeMapper;
using portculis::Intron;
using portculis::Junction;
using portculis::JunctionSystem;



namespace portculis {

const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;

typedef boost::error_info<struct PortculisError,string> PortculisErrorInfo;
struct PortculisException: virtual boost::exception, virtual std::exception { };

class Portculis {
private:

    // Can set these from the outside via the constructor
    string sortedBamFile;
    GenomeMapper* genomeMapper;
    string outputPrefix;
    uint16_t threads;
    bool verbose;
    
    // Sam header and refs info from the input bam
    SamHeader header;
    RefVector refs;
    
    // The set of distinct junctions found in the BAM file
    JunctionSystem junctionSystem;
    
    void init(  string _sortedBamFile, GenomeMapper* _genomeMapper, string _outputPrefix, 
                uint16_t _threads, bool _verbose) {
        
        sortedBamFile = _sortedBamFile;
        genomeMapper = _genomeMapper;
        outputPrefix = _outputPrefix;
        threads = _threads;
        verbose = _verbose;
        
        
        if (verbose)
            cerr << "Initialised Portculis instance" << endl;
    }
    

protected:

    string getAssociatedIndexFile(string bamFile) {
        return string(bamFile) + ".bti";
    }
       
    
    /**
     * Populates the set of distinct junctions.  
     * 
     * Also outputs all the unspliced alignments to a separate file if requested
     */
    void separateSplicedAlignments() {
        
        auto_cpu_timer timer(1, " = Wall time taken: %ws\n");
        
        BamReader reader;
        
        if (!reader.Open(sortedBamFile)) {
            BOOST_THROW_EXCEPTION(PortculisException() << PortculisErrorInfo(string(
                    "Could not open BAM reader for input: ") + sortedBamFile));
        }
        // Sam header and refs info from the input bam
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

       
        cout << " - Separating alignments from: " << sortedBamFile << endl;
        
        string indexFile = getAssociatedIndexFile(sortedBamFile);
        
        // Opens the index for this BAM file
        if ( !reader.OpenIndex(indexFile) ) {            
            BOOST_THROW_EXCEPTION(PortculisException() << PortculisErrorInfo(string(
                    "Could not open index for BAM: ") + indexFile));             
        }
        
        cout << " - Using BAM index: " << indexFile << endl;
        
        BamWriter unsplicedWriter;
        string unsplicedFile = getUnsplicedBamFile();

        if (!unsplicedWriter.Open(unsplicedFile, header, refs)) {
            BOOST_THROW_EXCEPTION(PortculisException() << PortculisErrorInfo(string(
                    "Could not open BAM writer for non-spliced file: ") + unsplicedFile));
        }

        cout << " - Saving unspliced alignments to: " << unsplicedFile << endl;
        
        BamAlignment al;
        uint64_t splicedCount = 0;
        uint64_t unsplicedCount = 0;
        uint64_t sumQueryLengths = 0;
        int32_t minQueryLength = 0;
        int32_t maxQueryLength = 100000;
        cout << " - Processing alignments ... ";
        cout.flush();
        while(reader.GetNextAlignment(al))
        {
            int32_t len = al.Length;
            minQueryLength = min(minQueryLength, len);
            maxQueryLength = max(maxQueryLength, len);
            
            sumQueryLengths += len;
            
            if (junctionSystem.addJunctions(al)) {
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
            BOOST_THROW_EXCEPTION(PortculisException() << PortculisErrorInfo(string(
                    "Could not open bam reader for unspliced alignments file: ") + unsplicedFile));
        }
        // Sam header and refs info from the input bam
        SamHeader header = reader.GetHeader();
        RefVector refs = reader.GetReferenceData();

        // Opens the index for this BAM file
        string unsplicedIndexFile = getAssociatedIndexFile(unsplicedFile);
        if ( !reader.OpenIndex(unsplicedIndexFile) ) {            
            if ( !reader.CreateIndex(BamIndex::BAMTOOLS) ) {
                BOOST_THROW_EXCEPTION(PortculisException() << PortculisErrorInfo(string(
                        "Error creating BAM index for unspliced alignments file: ") + unsplicedIndexFile));
            }            
        }
    }
    

public:

    Portculis() {
        init(   "",
                NULL, 
                DEFAULT_OUTPUT_PREFIX, 
                DEFAULT_THREADS,
                false
                );
    }
    Portculis(  string _sortedBam, GenomeMapper* _genomeMapper, 
                string _outputPrefix, uint16_t _threads, 
                bool _verbose) {
        init(  _sortedBam, _genomeMapper, _outputPrefix, 
                _threads, _verbose);
    }
    
    virtual ~Portculis() {        
    }
    
    
    string getUnsplicedBamFile() {
        return outputPrefix + ".unspliced.bam";
    }
    

    void process() {
       
        // Start the timer
        auto_cpu_timer timer(1, "Portculis finished.  Wall time taken: %ws\n");
        
        // Collect junctions from BAM file (also outputs unspliced alignments
        // to a separate file)
        cout << "Stage 1: Separating spliced alignments:" << endl;
        separateSplicedAlignments();        
        
        // Acquires donor / acceptor info from indexed genome file
        cout << "Stage 2: Scanning reference sequences:" << endl;
        uint64_t daSites = junctionSystem.scanReference(genomeMapper, refs);
        
        // Count the number of alignments found in upstream and downstream flanking 
        // regions for each junction
        cout << "Stage 3: Lookup unspliced alignments:" << endl;
        junctionSystem.findFlankingAlignments(getUnsplicedBamFile());
        
        cout << "Stage 4: Calculating junction status flags:" << endl;
        junctionSystem.calcJunctionStats();
        
        // Calculate all remaining metrics
        cout << "Stage 5: Calculating remaining junction metrics:" << endl;
        junctionSystem.calcAllRemainingMetrics();
        
        cout << "Stage 6: Outputting junction information:" << endl;
        junctionSystem.saveAll(outputPrefix);
    }
};
}
