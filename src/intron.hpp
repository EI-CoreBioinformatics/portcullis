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

#include <math.h>
#include <string>
#include <memory>

using std::min;
using std::string;
using std::shared_ptr;

#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
using boost::lexical_cast;

namespace portculis {    
    
enum Strand {
    POSITIVE,
    NEGATIVE,
    UNKNOWN
};
    
static Strand strandFromBool(bool reverseStrand) {
    return reverseStrand ? NEGATIVE : POSITIVE;
}

static Strand strandFromChar(char strand) {
    switch(strand) {
        case '+':
            return POSITIVE;
        case '-':
            return NEGATIVE;
        case '?':
            return UNKNOWN;
    }

    return UNKNOWN;
}

static char strandToChar(Strand strand) {
    
    switch(strand) {
        case POSITIVE:
            return '+';
        case NEGATIVE:
            return '-';
        case UNKNOWN:
            return '?';
    }

    return '?';
}

static string strandToString(Strand strand) {
    
    switch(strand) {
        case POSITIVE:
            return string("POSITIVE");
        case NEGATIVE:
            return string("NEGATIVE");
        case UNKNOWN:
            return string("UNKNOWN");
    }

    return string("UNKNOWN");
}
    
typedef boost::error_info<struct IntronError,string> IntronErrorInfo;
struct IntronException: virtual boost::exception, virtual std::exception { };


class Intron {
    
    
public:    
    int32_t refId;      // The reference sequence this location comes from
    int32_t start;      // The index of the base of the intron
    int32_t end;        // The index of the last base of the intron
    Strand strand;      // The strand the intron is on
    
    Intron() : Intron(-1, -1, -1, UNKNOWN) {}
    
    Intron(int32_t _refId, int32_t _start, int32_t _end, Strand _strand) :
        refId(_refId), start(_start), end(_end), strand(_strand) {
    }

    
    /**
     * Note that equality is determined purely on the ref id and the junction's
     * start and end... the flanking regions are ignored.
     * @param other
     * @return 
     */
    bool operator==(const Intron &other) const
    { 
        return (refId == other.refId &&
                start == other.start &&
                end == other.end &&
                strand == other.strand);
    }
    
    bool operator!=(const Intron &other) const
    { 
        return !((*this)==other);
    }
    
    
    int32_t size() const {
        return end - start + 1;
    }
    
    /**
     * We probably need to double check this logic.  We say this intron shares a
     * donor or acceptor with another intron, if the ref id is the same and we
     * find either the start or end position is the same. 
     *  
     * Note: We ignore the strand!  Is this right?
     * 
     * @param other
     * @return 
     */
    bool sharesDonorOrAcceptor(const Intron& other) {
        return refId == other.refId && 
                (start == other.start || end == other.end);
    }
    
    /**
     * Calculates the length of both anchors, based on the provided start and end
     * flanking regions and then returns the minimum of the two
     * @param leftAnchorStart The start position of the left anchor
     * @param rightAnchorEnd The end position of the right anchor (inclusive)
     * @return The minimum of the left anchor length and the right anchor length
     */
    int32_t minAnchorLength(int32_t leftAnchorStart, int32_t rightAnchorEnd) {
        
        if (leftAnchorStart >= start)
            BOOST_THROW_EXCEPTION(IntronException() << IntronErrorInfo(string(
                    "The intron start position must be greater than the left anchor start position")));
        
        if (rightAnchorEnd <= end)
            BOOST_THROW_EXCEPTION(IntronException() << IntronErrorInfo(string(
                    "The intron end position must be less than the right anchor end position")));
        
        int32_t lAnchor = start - leftAnchorStart;
        int32_t rAnchor = rightAnchorEnd - end;

        return min(lAnchor, rAnchor);        
    }
    
    void outputDescription(std::ostream &strm) {
        outputDescription(strm, "; ");
    }
    
    void outputDescription(std::ostream &strm, string delimiter) {
        strm << "RefId: " << refId << delimiter
             << "Start: " << start << delimiter
             << "End: " << end << delimiter
             << "Strand: " << strandToString(strand);
    }
    
    friend std::ostream& operator<<(std::ostream &strm, const Intron& l) {
        return strm << l.refId << "\t" << l.start << "\t" << l.end << "\t" << strandToChar(l.strand);
    }
    
    static string locationOutputHeader() {
        return string("refid\tstart\tend\tstrand"); 
    }
    
};

/**
    * Overload hash_value with Location so that boost know how to use this class
    * as a key in some kind of hash table
    * 
    * @param other
    * @return 
    */
struct IntronHasher {
    
    size_t operator()(const Intron& l) const
    {
       using boost::hash_value;
       using boost::hash_combine;

       // Start with a hash value of 0    .
       size_t seed = 0;

       // Modify 'seed' by XORing and bit-shifting in
       // one member of 'Key' after the other:
       hash_combine(seed, hash_value(l.refId));
       hash_combine(seed, hash_value(l.start));
       hash_combine(seed, hash_value(l.end));
       hash_combine(seed, hash_value(l.strand));

       // Return the result.
       return seed;
    }
};
}
