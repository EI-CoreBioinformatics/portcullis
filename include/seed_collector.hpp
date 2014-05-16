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

#include <iostream>
#include <vector>

#include <boost/foreach.hpp>

#include <api/BamReader.h>
#include <api/BamWriter.h>

using std::vector;

using namespace BamTools;

namespace portculis {

class CollectorResults {
private:
    void init(uint64_t _seedCount, uint64_t _unsplicedCount) {
        seedCount = _seedCount;
        unsplicedCount = _unsplicedCount;
    }

public:
    
    uint64_t seedCount;
    uint64_t unsplicedCount;
    
    CollectorResults() {
        init(0,0);
    }

    CollectorResults(uint64_t _seedCount, uint64_t _unsplicedCount) {
        init(_seedCount, _unsplicedCount);
    }

    void report(std::ostream& out) {
        out << endl 
            << "Seed collection report" << endl
            << "----------------------" << endl
            << "Found " << seedCount << " alignments containing one or more potential junctions" << endl
            << "Found " << unsplicedCount << " unspliced alignments" << endl << endl;
    }
};
    
class SeedCollector {
private:
 
    string sortedBam;    
    string seedFile;
    string unsplicedFile;
    bool verbose;
    
    BamReader reader;
    BamWriter seedWriter;
    BamWriter unsplicedWriter;
        
    SamHeader header;
    RefVector refs;
    
protected:
    
    bool isPotentialSeed(BamAlignment& al) {
        BOOST_FOREACH(CigarOp op, al.CigarData) {
            if (op.Type == 'N') {
                return true;
            }            
        }
        
        return false;        
    }
    
public:
    
    SeedCollector(string _sortedBam, string _seedFile, string _unsplicedFile, bool _verbose) : 
            sortedBam(_sortedBam), seedFile(_seedFile), unsplicedFile(_unsplicedFile), verbose(_verbose) {
        
        if (!reader.Open(sortedBam)) {
            throw "Could not open bam reader";
        }
        
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

       
        cout << "Loading alignments from: " << sortedBam << endl;
        
        if (verbose) {
            cout << "Header:" << endl << header.ToString() << endl;
            cout << "Refs:" << endl;

            BOOST_FOREACH(RefData ref, refs) {
               cout << ref.RefLength << " - " << ref.RefName << endl;
            }
        }      
        
        if (!seedWriter.Open(seedFile, header, refs)) {
            throw "Could not open bam writer for seed file";
        }
        
        if (!unsplicedWriter.Open(unsplicedFile, header, refs)) {
            throw "Could not open bam writer for unspliced file";
        }
        
        cout << "Sending seed alignments to: " << seedFile;
        cout << "Sending unspliced alignments to: " << unsplicedFile;
    }
            
    virtual ~SeedCollector() {
        reader.Close();
        seedWriter.Close();
        unsplicedWriter.Close();
    }
            
    void collect(CollectorResults& results) {
        
        BamAlignment al;
        uint64_t seedCount = 0;
        uint64_t unsplicedCount = 0;
        cout << "Processing alignments ... ";
        while(reader.GetNextAlignment(al))
        {
            if (isPotentialSeed(al)) {
                seedWriter.SaveAlignment(al);
                seedCount++;
            }
            else {
                unsplicedWriter.SaveAlignment(al);
                unsplicedCount++;
            }
        }
        cout << "done" << endl;
        
        results.seedCount = seedCount;
        results.unsplicedCount = unsplicedCount;
    }    
    
};
}

