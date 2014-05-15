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

#include <vector>

#include <boost/foreach.hpp>

#include <api/BamReader.h>
#include <api/BamWriter.h>

using std::vector;

using namespace BamTools;

namespace portculis {

class SeedCollector {
private:
 
    string sortedBam;    
    string seedFile;
    bool verbose;
    
    BamReader reader;
    BamWriter writer;
        
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
    
    SeedCollector(string _sortedBam, string _seedFile, bool _verbose) : 
            sortedBam(_sortedBam), seedFile(_seedFile), verbose(_verbose) {
        
        if (!reader.Open(sortedBam)) {
            throw "Could not open bam reader";
        }
        
        if (!writer.Open(seedFile, header, refs)) {
            throw "Could not open bam writer";
        }
        
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

        if (verbose) {
            
            cout << endl << "Collecting seeds from: " << sortedBam << endl;
            cout << "Header:" << endl << header.ToString() << endl;
            cout << "Refs:" << endl;

            BOOST_FOREACH(RefData ref, refs) {
               cout << ref.RefLength << " - " << ref.RefName << endl;
            }
        }
    }
            
    virtual ~SeedCollector() {
        reader.Close();
        writer.Close();
    }
            
    uint64_t collect() {
        
        BamAlignment al;
        uint64_t count = 0;
        while(reader.GetNextAlignment(al))
        {
            if (isPotentialSeed(al)) {
                writer.SaveAlignment(al);
                count++;
            }
        }
        
        return count;
    }
    
    
};
}

