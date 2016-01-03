//  ********************************************************************
//  This file is part of Portcullis.
//
//  Portcullis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portcullis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portcullis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

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

#include <portcullis/bam/bam_master.hpp>
using portcullis::bam::RefSeq;

#include <portcullis/intron.hpp>


portcullis::Intron::Intron(const Intron& other) {
    ref = other.ref;
    start = other.start;
    end = other.end;
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
bool portcullis::Intron::sharesDonorOrAcceptor(const Intron& other) {
    return ref.index == other.ref.index && 
            (start == other.start || end == other.end);
}

/**
 * Calculates the length of both anchors, based on the provided start and end
 * flanking regions and then returns the minimum of the two
 * @param leftAnchorStart The start position of the left anchor
 * @param rightAnchorEnd The end position of the right anchor (inclusive)
 * @return The minimum of the left anchor length and the right anchor length
 */
int32_t portcullis::Intron::minAnchorLength(int32_t leftAnchorStart, int32_t rightAnchorEnd) {

    if (leftAnchorStart > start)
        BOOST_THROW_EXCEPTION(IntronException() << IntronErrorInfo(string(
                "The intron start position must be greater than the left anchor start position: ") + 
                lexical_cast<string>(leftAnchorStart) + " **** " + 
                lexical_cast<string>(start) + "-" +
                lexical_cast<string>(end) + " **** " +
                lexical_cast<string>(rightAnchorEnd) + ") " +
                "Reference seq: " + this->ref.toString()));

    if (rightAnchorEnd < end)
        BOOST_THROW_EXCEPTION(IntronException() << IntronErrorInfo(string(
                "The intron end position must be less than the right anchor end position: (") + 
                lexical_cast<string>(leftAnchorStart) + " **** " + 
                lexical_cast<string>(start) + "-" +
                lexical_cast<string>(end) + " **** " +
                lexical_cast<string>(rightAnchorEnd) + ") " +
                "Reference seq: " + this->ref.toString()));

    int32_t lAnchor = start - leftAnchorStart;
    int32_t rAnchor = rightAnchorEnd - end;

    return min(lAnchor, rAnchor);        
}


void portcullis::Intron::outputDescription(std::ostream &strm, string delimiter) {
    strm << "RefId: " << ref.index << delimiter
         << "RefName: " << ref.name << delimiter
         << "RefLength: " << ref.length << delimiter
         << "Start: " << start << delimiter
         << "End: " << end;
}

size_t portcullis::IntronHasher::operator()(const Intron& l) const
{
   using boost::hash_value;
   using boost::hash_combine;

   // Start with a hash value of 0    .
   size_t seed = 0;

   // Modify 'seed' by XORing and bit-shifting in
   // one member of 'Key' after the other:
   hash_combine(seed, hash_value(l.ref.index));
   hash_combine(seed, hash_value(l.start));
   hash_combine(seed, hash_value(l.end));

   // Return the result.
   return seed;
}

bool portcullis::IntronComparator::operator()(const Intron& l, const Intron& r) const
{
    if (l.ref.index < r.ref.index) {
        return true;
    }
    else if (l.ref.index == r.ref.index) {
        if (l.start < r.start) {
            return true;
        }
        else if (l.start == r.start) {
            if (l.end < r.end) {
                return true;
            }        
        }
    }
    
    return false;
}

