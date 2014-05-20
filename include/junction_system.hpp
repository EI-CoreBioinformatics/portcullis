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

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "location.hpp"
#include "junction.hpp"
#include "genome_mapper.hpp"

using boost::lexical_cast;
using boost::shared_ptr;

using portculis::Location;
using portculis::Junction;

typedef boost::unordered_map<Location, shared_ptr<Junction> > DistinctJunctions;
typedef boost::unordered_map<Location, shared_ptr<Junction> >::iterator JunctionMapIterator;
typedef std::pair<const Location, shared_ptr<Junction> > JunctionMapType;
typedef std::vector<shared_ptr<Junction> > JunctionList;

namespace portculis {
class JunctionSystem {

private:    
    DistinctJunctions distinctJunctions;
    JunctionList junctionList;
    
public:
    
    
    JunctionSystem() {        
    }
    
    virtual ~JunctionSystem() {
             
    }
    
    
    size_t size() {
        return distinctJunctions.size();
    }
    
    /**
     * Adds any new junctions found from the given alignment to the set managed 
     * by this class
     * @param al The alignment to search for junctions
     * @return Whether a junction was found in this alignment or not
     */
    bool addJunctions(const BamAlignment& al) {
        return addJunctions(al, 0, 0);
    }
    
    bool addJunctions(const BamAlignment& al, const size_t startOp, const int32_t offset) {
        
        bool foundJunction = false;
        
        size_t nbOps = al.CigarData.size();
        
        int32_t refId = al.RefID;
        int32_t lStart = al.Position + offset;        
        int32_t lEnd = lStart;
        int32_t rStart = lStart;
        int32_t rEnd = lStart;
        
        for(size_t i = startOp; i < nbOps; i++) {
            
            CigarOp op = al.CigarData[i];
            if (op.Type == 'N') {
                foundJunction = true;
                
                rStart = lEnd + op.Length;
                rEnd = rStart;
                
                size_t j = i+1;
                while (j < nbOps && al.CigarData[j].Type != 'N') {
                    rEnd += al.CigarData[j++].Length;
                }
                
                Location* location = new Location(refId, lStart, lEnd, rStart, rEnd);
                
                // We should now have the complete junction location information
                JunctionMapIterator it = distinctJunctions.find(*location);
                
                // If we couldn't find this location in the hashmap, add a new
                // location / junction pair.  If we've seen this location before
                // then add this alignment to the existing junction
                if (it == distinctJunctions.end()) {
                    shared_ptr<Junction> junction(new Junction(location));
                    junction->addJunctionAlignment(al);
                    distinctJunctions[*location] = junction;
                    junctionList.push_back(junction);
                }
                else {
                    it->second->addJunctionAlignment(al);
                    delete location;    // We don't need this location as we already have one
                }
                
                // Check if we have fully processed the cigar or not.  If not, then
                // that means that this cigar contains additional junctions, so 
                // process those using recursion
                if (j < nbOps) {
                    
                    addJunctions(al, i+1, rStart);
                    break;
                }                
            }
            else {
                lEnd += op.Length;                
            }            
        }
        
        return foundJunction;        
    }
    
    void useGenome(GenomeMapper* genomeMapper, RefVector& refs) {
        
        BOOST_FOREACH(shared_ptr<Junction> j, junctionList) {
            
            Location* loc = j->getLocation();
            
            if (loc != NULL) {
                int32_t refid = loc->refId;
                int32_t intronStart = loc->lEnd;
                int32_t intronEnd = loc->rStart - 1;
                
                char* refName = new char[refs[refid].RefName.size() + 1];
                strcpy(refName, refs[refid].RefName.c_str());
                int seq1Len = -1;
                int seq2Len = -1;

                char* seq1 = genomeMapper->fetch(refName, intronStart, intronStart+1, &seq1Len);
                char* seq2 = genomeMapper->fetch(refName, intronEnd-1, intronEnd, &seq2Len);
                
                bool hasDA = (seq1Len == 2 && seq2Len == 2) ? 
                    hasDonorAcceptor(seq1, seq2) :
                    false;
                
                j->setDonorAndAcceptorMotif(hasDA);
            }
        }
    }
    
    bool hasDonorAcceptor(char seq1[], char seq2[]) {
        return ((seq1[0] == 'G' && seq1[1] == 'T' && seq2[0] == 'A' && seq2[1] == 'G') ||
                (seq1[0] == 'C' && seq1[1] == 'A' && seq2[0] == 'T' && seq2[1] == 'C'));
    }
    
    void outputDescription(std::ostream &strm) {
        
        uint64_t i = 0;
        BOOST_FOREACH(shared_ptr<Junction> j, junctionList) {
            strm << "Junction " << i++ << ":" << endl;
            j->outputDescription(strm);
            strm << endl;
        }        
    }
    
    friend std::ostream& operator<<(std::ostream &strm, const JunctionSystem& js) {
        
        strm << "index\t" << Junction::junctionOutputHeader() << endl;
        
        uint64_t i = 0;
        BOOST_FOREACH(shared_ptr<Junction> j, js.junctionList) {
            strm << i++ << "\t" << *j << endl;
        }
        
        return strm;
    }
    
};
}
