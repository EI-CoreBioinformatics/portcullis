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
#include <tr1/unordered_set>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>

#include <api/BamReader.h>

#include "junction.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::tr1::unordered_set;
using std::vector;

using boost::lexical_cast;

using portculis::junction::Junction;
//using portculis::junction::JunctionHasher;
using portculis::junction::Location;

using namespace BamTools;

namespace portculis {

const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;


class Portculis {
private:

    // Can set these from the outside via the constructor
    string seedFile;
    string sortedBamFile;
    string genomeFile;
    string outputPrefix;
    uint16_t threads;
    bool verbose;
    
    //unordered_set<Location,LocationHasher> junctionLocations;
    
    void init(  string _seedFile, string _sortedBamFile, string _genomeFile, 
                string _outputPrefix, uint16_t _threads, bool _verbose) {
        
        seedFile = _seedFile;
        sortedBamFile = _sortedBamFile;
        genomeFile = _genomeFile;
        outputPrefix = _outputPrefix;
        threads = _threads;
        verbose = _verbose;
        
        if (verbose)
            cerr << "Initialised Portculis instance" << endl;
    }
    

protected:

    
    void calcJunction(BamAlignment& al, Location loc) {
        
        loc.refId = al.RefID;
        loc.lStart = al.Position;
        
        int32_t lEnd = loc.lStart;
        int32_t rEnd = loc.lStart;
        bool lEndFound = false;
        bool rEndFound = false;
        BOOST_FOREACH(CigarOp op, al.CigarData) {
            
            if (!lEndFound && op.Type != 'N') {
                lEnd += op.Length;                
            }
            else if (!lEndFound && op.Type == 'N') {
                lEndFound = true;
                loc.lEnd = lEnd;
                loc.rStart = lEnd + op.Length;
                rEnd = loc.rStart;
            }
            else if (lEndFound && !rEndFound && op.Type != 'N') {
                rEnd += op.Length;
            }
            else if (lEndFound && !rEndFound && op.Type == 'N') {
                rEndFound = true;
                loc.rEnd = rEnd;
                //loc.doubleJunction = true;
            }
        }
        
        if (!rEndFound) {
            loc.rEnd = rEnd;
        }
        
        //junctions[loc];
    }

public:

    Portculis() {
        init(   "",
                "", 
                "", 
                DEFAULT_OUTPUT_PREFIX, 
                DEFAULT_THREADS, 
                false
                );
    }
    Portculis(  string _seedFile, string _sortedBam, string _genomeFile, 
                string _outputPrefix, uint16_t _threads, bool _verbose) {
        init(   _seedFile, _sortedBam, _genomeFile, 
                _outputPrefix, _threads, _verbose);
    }
    
    virtual ~Portculis() {
        /*if (splicedOut != NULL)
            delete splicedOut;

        if (weirdOut != NULL)
            delete weirdOut;*/
    }

    void process() {
       
        BamReader reader;
        if (!reader.Open(seedFile)) {
            throw "Could not open bam reader for seed file";
        }
        
        BamAlignment al;
        while(reader.GetNextAlignment(al)) {
            
            junction::Location loc;
            calcJunction(al, loc);
            cout << loc.toString(true) << endl;            
        }
        
        reader.Close();        
    }
};
}
