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

#include <boost/lexical_cast.hpp>

using std::string;

using boost::lexical_cast;

namespace portculis {    

class Location {
public:    
    int32_t refId;
    int32_t lStart;
    int32_t lEnd;
    int32_t rStart;
    int32_t rEnd;
    
    Location() :
        refId(-1), lStart(-1), lEnd(-1), rStart(-1), rEnd(-1) {        
    }
    
    Location(int32_t _refId, int32_t _lStart, int32_t _lEnd, int32_t _rStart, int32_t _rEnd) :
        refId(_refId), lStart(_lStart), lEnd(_lEnd), rStart(_rStart), rEnd(_rEnd) {
    }

    
    friend size_t hash_value(const Location& l)
    {
        // Start with a hash value of 0    .
        std::size_t seed = 0;

        // Modify 'seed' by XORing and bit-shifting in
        // one member of 'Key' after the other:
        boost::hash_combine(seed, l.refId);
        boost::hash_combine(seed, l.lEnd);
        boost::hash_combine(seed, l.rStart);

        // Return the result.
        return seed;
    }
    
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
}
