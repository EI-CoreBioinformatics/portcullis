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
    bool outputUnspliced;
    bool verbose;
    
    BamReader reader;
    BamWriter unsplicedWriter;
     
    // Sam header and refs info from the input bam
    SamHeader header;
    RefVector refs;
    
    // The set of distinct junctions found in the BAM file
    JunctionSystem junctionSystem;
    
    void init(  string _sortedBamFile, GenomeMapper* _genomeMapper, string _outputPrefix, 
                uint16_t _threads, bool _outputUnspliced, bool _verbose) {
        
        sortedBamFile = _sortedBamFile;
        genomeMapper = _genomeMapper;
        outputPrefix = _outputPrefix;
        threads = _threads;
        outputUnspliced = _outputUnspliced;
        verbose = _verbose;
        
        if (!reader.Open(sortedBamFile)) {
            throw "Could not open bam reader";
        }
        
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

       
        cout << "Will load alignments from: " << sortedBamFile << endl;
        
        if (verbose) {
            cout << "Header:" << endl << header.ToString() << endl;
            cout << "Refs:" << endl;

            BOOST_FOREACH(RefData ref, refs) {
               cout << ref.RefLength << " - " << ref.RefName << endl;
            }
        }
        
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
        
        if (outputUnspliced) {
            string unsplicedFile = getUnsplicedFile();
        
            if (!unsplicedWriter.Open(unsplicedFile, header, refs)) {
                throw "Could not open bam writer for unspliced file";
            }
            
            cout << "Sending unspliced alignments to: " << unsplicedFile << endl;
        }
            
        BamAlignment al;
        uint64_t splicedCount = 0;
        uint64_t unsplicedCount = 0;
        cout << "Processing alignments ... ";
        cout.flush();
        while(reader.GetNextAlignment(al))
        {
            if (junctionSystem.addJunctions(al)) {
                splicedCount++;
            }
            else if (outputUnspliced) {
                unsplicedWriter.SaveAlignment(al);
                unsplicedCount++;
            }
        }
        cout << "done" << endl;
        
        if (outputUnspliced) {
            
            cout << "Found " << junctionSystem.size() << " junctions from " << splicedCount << " spliced alignments." << endl;
            cout << "Found " << unsplicedCount << " unspliced alignments." << endl;
            unsplicedWriter.Close();
        }        
    }
    

public:

    Portculis() {
        init(   "",
                NULL, 
                DEFAULT_OUTPUT_PREFIX, 
                DEFAULT_THREADS,
                false,
                false
                );
    }
    Portculis(  string _sortedBam, GenomeMapper* _genomeMapper, 
                string _outputPrefix, uint16_t _threads, 
                bool _outputUnspliced, bool _verbose) {
        init(  _sortedBam, _genomeMapper, _outputPrefix, 
                _threads, _outputUnspliced, _verbose);
    }
    
    virtual ~Portculis() {
        
        reader.Close();
        
        
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
        uint64_t daSites = junctionSystem.findDonorAcceptorSites(genomeMapper, refs);
        cout << "done" << endl
             << "Found " << daSites << " valid donor / acceptor sites." << endl;
                
        // Calculate all remaining metrics
        junctionSystem.calcAllMetrics();
        
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
