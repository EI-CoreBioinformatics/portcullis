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

using std::endl;
using std::string;
using std::size_t;
using std::vector;

using boost::lexical_cast;

using namespace BamTools;

using portculis::Location;

namespace portculis {    

    
class Junction {
private:
    
    // **** Properties that describe where the junction is ****
    Location* location;
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
    
    Junction() :
        location(NULL) {
        init();
    }
    
    Junction(Location* _location) :
        location(_location) {
        init();
    }
    
    // **** Destructor ****
    virtual ~Junction() {
        
    }
    
   
    Location* getLocation() const {
        return location;
    }

    void setLocation(Location* location) {
        this->location = location;
    }

    
    
    void addJunctionAlignment(const BamAlignment& al) {
        this->junctionAlignments.push_back(al);                        
    }
    
    void setDonorAndAcceptorMotif(bool donorAndAcceptorMotif) {
        this->donorAndAcceptorMotif = donorAndAcceptorMotif;
    }

    
    
    // **** Metric getters ****
    
    uint16_t getMaxMinAnchor() const {
        return maxMinAnchor;
    }

    /**
     * Metric 1: The number of alignments directly supporting this junction
     * @return 
     */
    size_t getNbJunctionAlignments() const {
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
        return location != NULL ? location->rStart - location->lEnd : 0;
    }
    
    /**
     * Metric 10: The number of upstream non-spliced supporting reads
     * @return 
     */
    size_t getNbUpstream() const {
        return leftFlankingAlignments.size();
    }
    
    /**
     * Metric 10: The number of upstream non-spliced supporting reads
     * @return 
     */
    size_t getNbDownstream() const {
        return rightFlankingAlignments.size();
    }
    
    
    
    void outputDescription(std::ostream &strm) {
        strm << "Location: ";
        
        if (location != NULL) {
            location->outputDescription(strm);
        }
        else {
            strm << "No location set";
        }
        
        strm << endl
             << "1:  # Junction Alignments: " << getNbJunctionAlignments() << endl
             << "2:  Has Donor + Acceptor Motif: " << donorAndAcceptorMotif << endl
             << "3:  Intron Size: " << getIntronSize() << endl
             << "10: # Upstream Non-Spliced Alignments: " << getNbUpstream() << endl
             << "11: # Downstream Non-Spliced Alignments: " << getNbDownstream() << endl;
                
    }
    
    friend std::ostream& operator<<(std::ostream &strm, const Junction& j) {
        return strm << *(j.location) << "\t" 
                    << j.getNbJunctionAlignments() << "\t"
                    << j.donorAndAcceptorMotif << "\t"
                    << j.getIntronSize() << "\t"
                    << j.getNbUpstream() << "\t"
                    << j.getNbDownstream();
    }
    
    static string junctionOutputHeader() {
        return string(Location::locationOutputHeader()) + string("\tM1\tM2\tM3\tM10\tM11"); 
    }

};

}
