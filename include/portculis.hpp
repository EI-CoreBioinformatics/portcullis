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

#include "location.hpp"
#include "junction.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

using boost::lexical_cast;

using portculis::Location;
using portculis::Junction;

using namespace BamTools;

typedef boost::unordered_map<Location, Junction> DistinctJunctions;
typedef boost::unordered_map<Location, Junction>::iterator JunctionMapIterator;
typedef std::vector<Junction*> JunctionList;

namespace portculis {

const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;


class Portculis {
private:

    // Can set these from the outside via the constructor
    string sortedBamFile;
    string genomeFile;
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
    DistinctJunctions distinctJunctions;
    JunctionList junctionList;
    
    void init(  string _sortedBamFile, string _genomeFile, string _outputPrefix, 
                uint16_t _threads, bool _outputUnspliced, bool _verbose) {
        
        sortedBamFile = _sortedBamFile;
        genomeFile = _genomeFile;
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
     * Adds any new junctions found from the given alignment to the set managed 
     * by this class
     * @param al The alignment to search for junctions
     * @return Whether a junction was found in this alignment or not
     */
    bool createJunctions(const BamAlignment& al) {
        return createJunctions(al, 0, 0);
    }
    
    bool createJunctions(const BamAlignment& al, const size_t startOp, const int32_t offset) {
        
        bool foundJunction = false;
        
        size_t nbOps = al.CigarData.size();
        
        int32_t refId = al.RefID;
        int32_t lStart = al.Position + offset;        
        int32_t lEnd = lStart;
        int32_t rStart = lStart;
        int32_t rEnd = lStart;
        
        for(size_t i = startOp; i < nbOps; i++) {
            
            CigarOp op = al.CigarData[i];
            if (op.Type == 'N') {
                foundJunction = true;
                
                rStart = lEnd + op.Length;
                rEnd = rStart;
                
                size_t j = i+1;
                while (j < nbOps && al.CigarData[j].Type != 'N') {
                    rEnd += al.CigarData[j++].Length;
                }
                
                Location location(refId, lStart, lEnd, rStart, rEnd);
                Junction junction(location);
                
                // We should now have the complete junction location information
                JunctionMapIterator it = distinctJunctions.find(location);
                
                if (it == distinctJunctions.end()) {
                    distinctJunctions[location] = junction;
                    junctionList.push_back(&junction);
                }
                else {
                    it->second.addJunctionAlignment(al);
                }
                
                // Check if we have fully processed the cigar or not.  If not, then
                // that means that this cigar contains additional junctions, so 
                // process those using recursion
                if (j < nbOps) {
                    
                    createJunctions(al, i+1, rStart);
                    break;
                }                
            }
            else {
                lEnd += op.Length;                
            }            
        }
        
        return foundJunction;        
    }
    
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
            
            cout << "Sending unspliced alignments to: " << unsplicedFile;
        }
            
        BamAlignment al;
        uint64_t splicedCount = 0;
        uint64_t unsplicedCount = 0;
        cout << "Processing alignments ... ";
        while(reader.GetNextAlignment(al))
        {
            if (createJunctions(al)) {
                splicedCount++;
            }
            else if (outputUnspliced) {
                unsplicedWriter.SaveAlignment(al);
                unsplicedCount++;
            }
        }
        cout << "done" << endl;
        
        if (outputUnspliced) {
            
            cout << "Found " << distinctJunctions.size() << " junctions from " << splicedCount << " spliced alignments." << endl;
            cout << "Found " << unsplicedCount << " unspliced alignments." << endl;
            unsplicedWriter.Close();
        }        
    }    
    

public:

    Portculis() {
        init(   "",
                "", 
                DEFAULT_OUTPUT_PREFIX, 
                DEFAULT_THREADS,
                false,
                false
                );
    }
    Portculis(  string _sortedBam, string _genomeFile, 
                string _outputPrefix, uint16_t _threads, 
                bool _outputUnspliced, bool _verbose) {
        init(  _sortedBam, _genomeFile, _outputPrefix, 
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
    }
};
}
