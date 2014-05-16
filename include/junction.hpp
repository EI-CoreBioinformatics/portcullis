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

#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>

using std::string;
using std::size_t;

using boost::lexical_cast;

namespace portculis {    
namespace junction {

/*
struct LocationHasher
{
  size_t operator()(const Location& l) const
  {
      using boost::hash_value;
      using boost::hash_combine;

      // Start with a hash value of 0    .
      std::size_t seed = 0;

      // Modify 'seed' by XORing and bit-shifting in
      // one member of 'Key' after the other:
      hash_combine(seed,hash_value(l.refId));
      hash_combine(seed,hash_value(l.lEnd));
      hash_combine(seed,hash_value(l.rStart));
      
      // Return the result.
      return seed;
  }
};
    
struct JunctionHasher
{
  size_t operator()(const Junction& j) const
  {
      using boost::hash_value;
      using boost::hash_combine;

      // Start with a hash value of 0    .
      std::size_t seed = 0;

      // Modify 'seed' by XORing and bit-shifting in
      // one member of 'Key' after the other:
      hash_combine(seed,hash_value(j.location));
      
      // Return the result.
      return seed;
  }
};*/
    
class Location {
public:    
    int32_t refId;
    int32_t lStart;
    int32_t lEnd;
    int32_t rStart;
    int32_t rEnd;
    
    /**
     * Note that equality is determined purely on the ref id and the junction's
     * start and end... the flanking regions are ignored.
     * @param other
     * @return 
     */
    bool operator==(const Location &other) const
    { 
        return (refId == other.refId &&
                lEnd == other.lEnd &&
                rStart == other.rStart);
    }
    
    string toString() {
        return  string("RefId: ") + lexical_cast<string>(refId) + 
                string("; lStart: ") + lexical_cast<string>(lStart) + 
                string("; lEnd: ") + lexical_cast<string>(lEnd) + 
                string("; rStart: ") + lexical_cast<string>(rStart) + 
                string("; rEnd: ") + lexical_cast<string>(rEnd);
    }
    
    string toString(bool tabSeparated) {
        return tabSeparated ? 
            (
                lexical_cast<string>(refId) + string("\t") + 
                lexical_cast<string>(lStart) + string("\t") + 
                lexical_cast<string>(lEnd) + string("\t") + 
                lexical_cast<string>(rStart) + string("\t") + 
                lexical_cast<string>(rEnd) + string("\t")
            ) 
            : toString();
    }
};
    
class Junction {
private:
    
    // **** Properties that describe where the junction is ****
    Location location;
    
    // **** Junction metrics ****
    uint32_t nbReads;                   // Metric 1
    bool hasDonorAndAcceptorMotif;      // Metric 2
                                        // Metric 3 (intron size) calculated via location properties
    uint16_t maxMinAnchor;              // Metric 4
    
    
    
    void init(Location _location) {
        location = _location;
        
        nbReads = 0;
        hasDonorAndAcceptorMotif = false;
        maxMinAnchor = 0;
        
    }
    
public:
    
    // **** Constructors ****
    
    Junction() {
        init(Location());
    }
    
    Junction(Location location) {
        init(location);
    }
    
    /**
     * This is purely based on the location, not the metrics.  This is to enable
     * us to handle junctions uniquely in a hash table / set
     * @param other
     * @return 
     */
    bool operator==(const Junction &other) const
    { 
        return location == other.location;
    }
    
    // **** Metric getters ****
    
    bool isHasDonorAndAcceptorMotif() const {
        return hasDonorAndAcceptorMotif;
    }

    Location getLocation() const {
        return location;
    }

    uint16_t getMaxMinAnchor() const {
        return maxMinAnchor;
    }

    uint32_t getNbReads() const {
        return nbReads;
    }
    
    int32_t getIntronSize() {
        return location.rStart - location.lEnd;
    }
    
    
    
    string toString() {
        return  location.toString();
    }
    
    string toString(bool tabSeparated) {
        return tabSeparated ? 
            (
                location.toString(true)
            ) 
            : toString();
    }
};
}
}