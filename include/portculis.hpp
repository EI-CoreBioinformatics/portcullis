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

#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "genome_mapper.hpp"
#include "location.hpp"
#include "junction.hpp"
#include "junction_system.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

using boost::lexical_cast;

using portculis::GenomeMapper;
using portculis::Location;
using portculis::Junction;
using portculis::JunctionSystem;

using namespace BamTools;

namespace portculis {

const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;


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

    
    /**
     * Populates the set of distinct junctions.  
     * 
     * Also outputs all the unspliced alignments to a separate file if requested
     */
    void collect() {
        
        BamReader reader;
        
        if (!reader.Open(sortedBamFile)) {
            throw "Could not open bam reader";
        }
        // Sam header and refs info from the input bam
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

       
        cout << "Will load alignments from: " << sortedBamFile << endl;
        
        string indexFile = sortedBamFile + ".bti";
        
        // Opens the index for this BAM file
        if ( !reader.OpenIndex(indexFile) ) {
            
            cerr << "WARNING: Couldn't find suitable index file for: " << sortedBamFile << endl
                 << "         This should not have happened.  Creating index ...";
            cerr.flush();
        
            if ( !reader.CreateIndex(BamIndex::BAMTOOLS) ) {
               throw "Could not create BAM index"; 
            }
            cerr << "done." << endl;
            cerr.flush();
        }
        
        BamWriter unsplicedWriter;
        string unsplicedFile = getUnsplicedFile();

        if (!unsplicedWriter.Open(unsplicedFile, header, refs)) {
            throw "Could not open bam writer for unspliced file";
        }

        cout << "Sending unspliced alignments to: " << unsplicedFile << endl;
            
        BamAlignment al;
        uint64_t splicedCount = 0;
        uint64_t unsplicedCount = 0;
        uint64_t sumQueryLengths = 0;
        cout << "Processing alignments ... ";
        cout.flush();
        while(reader.GetNextAlignment(al))
        {
            sumQueryLengths += al.Length;
            
            if (junctionSystem.addJunctions(al)) {
                splicedCount++;
            }
            else {
                unsplicedWriter.SaveAlignment(al);
                unsplicedCount++;
            }
        }
        cout << "done." << endl;
        
        double meanQueryLength = (double)sumQueryLengths / (double)(splicedCount + unsplicedCount);
        junctionSystem.setMeanQueryLength(meanQueryLength);
        
        cout << "Found " << junctionSystem.size() << " junctions from " << splicedCount << " spliced alignments." << endl;
        cout << "Found " << unsplicedCount << " unspliced alignments." << endl;
        unsplicedWriter.Close();
        
        // Reset the reader in case anyone else want to use it later
        reader.Close();
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
    
    
    string getUnsplicedFile() {
        return string(outputPrefix) + ".unspliced.bam";
    }

    void process() {
       
        // Collect junctions from BAM file (also outputs unspliced alignments
        // to a separate file)
        collect();
        
        // Acquires donor / acceptor info from indexed genome file
        cout << "Acquiring donor / acceptor sites from genome ... ";
        cout.flush();
        uint64_t daSites = junctionSystem.findDonorAcceptorSites(genomeMapper, refs);
        cout << "done." << endl
             << "Found " << daSites << " valid donor / acceptor sites." << endl;
        
        // Count the number of alignments found in upstream and downstream flanking 
        // regions for each junction
        cout << "Acquiring non-spliced alignments from flanking windows ... ";
        cout.flush();
        junctionSystem.findFlankingAlignments(getUnsplicedFile());
        cout << "done." << endl;
        
        
        // Calculate all remaining metrics
        cout << "Calculating remaining junction metrics ... ";
        cout.flush();
        junctionSystem.calcAllMetrics();
        cout << "done." << endl;
        
        
        // Print descriptive output to file
        string junctionReportPath = outputPrefix + ".junctions.txt";
        ofstream junctionReportStream;
        junctionReportStream.open (junctionReportPath.c_str());
        junctionSystem.outputDescription(junctionReportStream);
        junctionReportStream.close();

        // Print junction stats to file
        string junctionFilePath = outputPrefix + ".junctions.tab";
        ofstream junctionFileStream;
        junctionFileStream.open (junctionFilePath.c_str());
        junctionFileStream << junctionSystem << endl;
        junctionFileStream.close();
        
    }
};
}
