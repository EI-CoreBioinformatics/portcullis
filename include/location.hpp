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
using std::string;

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
using boost::lexical_cast;
using boost::shared_ptr;

namespace portculis {    
    
enum Strand {
    POSITIVE,
    NEGATIVE
};
    
static Strand strandFromBool(bool reverseStrand) {
    return reverseStrand ? NEGATIVE : POSITIVE;
}

static char strandToChar(Strand strand) {
    return strand == POSITIVE ? '+' : '-';
}

static string strandToString(Strand strand) {
    return strand == POSITIVE ? string("POSITIVE") : string("NEGATIVE");
}
    
class Location {
    
    
public:    
    int32_t refId;
    int32_t start;
    int32_t end;
    Strand strand;
    
    Location() :
        refId(-1), start(-1), end(-1), strand(POSITIVE) {        
    }
    
    Location(int32_t _refId, int32_t _start, int32_t _end, Strand _strand) :
        refId(_refId), start(_start), end(_end), strand(_strand) {
    }

    
    friend size_t hash_value(const Location& l)
    {
        // Start with a hash value of 0    .
        std::size_t seed = 0;

        // Modify 'seed' by XORing and bit-shifting in
        // one member of 'Key' after the other:
        boost::hash_combine(seed, l.refId);
        boost::hash_combine(seed, l.start);
        boost::hash_combine(seed, l.end);
        boost::hash_combine(seed, l.strand);

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
                start == other.start &&
                end == other.end &&
                strand == other.strand);
    }
    
    bool operator!=(const Location &other) const
    { 
        return !((*this)==other);
    }
    
    
    int32_t size() const {
        return end - start + 1;
    }
    
    
    void outputDescription(std::ostream &strm) {
        strm << "RefId: " << refId
             << "; Start: " << start
             << "; End: " << end
             << "; Strand: " << strandToString(strand);
    }
    
    friend std::ostream& operator<<(std::ostream &strm, const Location& l) {
        return strm << l.refId << "\t" << l.start << "\t" << l.end << "\t" << strandToChar(l.strand);
    }
    
    static string locationOutputHeader() {
        return string("refid\tstart\tend\tstrand"); 
    }
    
};

}
