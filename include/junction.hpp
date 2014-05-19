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

#include <string>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>

#include <api/BamAlignment.h>

#include "location.hpp"

using std::string;
using std::size_t;
using std::vector;

using boost::lexical_cast;

using namespace BamTools;

namespace portculis {    

    
class Junction {
private:
    
    // **** Properties that describe where the junction is ****
    Location location;
    vector<BamAlignment> junctionAlignments;
    vector<BamAlignment> leftFlankingAlignments;
    vector<BamAlignment> rightFlankingAlignments;
    
    
    // **** Junction metrics ****
                                        // Metric 1 (nbReads) derived from size of junction alignment vector
    bool donorAndAcceptorMotif;         // Metric 2 calculated from 
                                        // Metric 3 (intron size) calculated via location properties
    uint16_t maxMinAnchor;              // Metric 4
    
    
    
    void init() {
        
        donorAndAcceptorMotif = false;
        maxMinAnchor = 0;
    }
    
public:
    
    // **** Constructors ****
    
    Junction() {
        init();
    }
    
    Junction(Location _location) :
        location(_location) {
        init();
    }
    
    // **** Destructor ****
    virtual ~Junction() {
        
    }
    
   
    Location getLocation() const {
        return location;
    }

    void setLocation(Location location) {
        this->location = location;
    }

    
    
    void addJunctionAlignment(const BamAlignment& al) {
        this->junctionAlignments.push_back(al);                        
    }
    
    
    // **** Metric getters ****
    
    uint16_t getMaxMinAnchor() const {
        return maxMinAnchor;
    }

    /**
     * Metric 1: The number of alignments directly supporting this junction
     * @return 
     */
    uint32_t getNbJunctionAlignments() const {
        return junctionAlignments.size();
    }
    
    /**
     * Metric 2: Whether or not there is a donor and acceptor motif at the two base
     * pairs at the start and end of the junction / intron
     * @return 
     */
    bool hasDonorAndAcceptorMotif() const {
        return donorAndAcceptorMotif;
    }
    
    /**
     * Metric 3: The intron size
     * @return 
     */
    int32_t getIntronSize() const {
        return location.rStart - location.lEnd;
    }
    
    
    
    string toString() const {
        return  string("RefId: ") + lexical_cast<string>(location.refId) + 
                string("; lStart: ") + lexical_cast<string>(location.lStart) + 
                string("; lEnd: ") + lexical_cast<string>(location.lEnd) + 
                string("; rStart: ") + lexical_cast<string>(location.rStart) + 
                string("; rEnd: ") + lexical_cast<string>(location.rEnd);
    }
    
    string toString(bool tabSeparated) const {
        return tabSeparated ? 
            (
                lexical_cast<string>(location.refId) + string("\t") + 
                lexical_cast<string>(location.lStart) + string("\t") + 
                lexical_cast<string>(location.lEnd) + string("\t") + 
                lexical_cast<string>(location.rStart) + string("\t") + 
                lexical_cast<string>(location.rEnd) + string("\t")
            ) 
            : toString();
    }
};
}
