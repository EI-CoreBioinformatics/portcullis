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

#pragma once

#include <fstream>
#include <vector>
#include <memory>
#include <unordered_map>
using std::ofstream;
using std::shared_ptr;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

#include "intron.hpp"
#include "junction.hpp"
#include "seq_utils.hpp"
using portcullis::Intron;
using portcullis::IntronHasher;
using portcullis::Junction;
using portcullis::JunctionPtr;
using portcullis::SeqUtils;

typedef std::unordered_map<Intron, JunctionPtr, IntronHasher> DistinctJunctions;
typedef std::unordered_map<Intron, JunctionPtr, IntronHasher>::iterator JunctionMapIterator;
typedef std::pair<const Intron, JunctionPtr> JunctionMapType;
typedef std::vector<JunctionPtr> JunctionList;
typedef std::shared_ptr<JunctionList> JunctionListPtr;

namespace portcullis {
class JunctionSystem {

private:    
    DistinctJunctions distinctJunctions;
    JunctionList junctionList;
    
    int32_t minQueryLength;
    double meanQueryLength;
    int32_t maxQueryLength;

    vector<RefSeq> refs;
    
    size_t createJunctionGroup(size_t index, vector<JunctionPtr>& group);
   
    void findJunctions(const int32_t refId, JunctionList& subset);
    
    
public:
    
    static string version;
    
    JunctionSystem();
    
    JunctionSystem(vector<RefSeq> refs);
    
    JunctionSystem(path junctionFile);
    
    virtual ~JunctionSystem();
    
    const JunctionList& getJunctions() const;
    
    size_t size();
    
    double getMeanQueryLength() const {
        return meanQueryLength;
    }

    void setMeanQueryLength(double meanQueryLength) {
        this->meanQueryLength = meanQueryLength;
    }
    
    int32_t getMaxQueryLength() const {
        return maxQueryLength;
    }

    void setMaxQueryLength(int32_t maxQueryLength) {
        this->maxQueryLength = maxQueryLength;
    }

    int32_t getMinQueryLength() const {
        return minQueryLength;
    }

    void setMinQueryLength(int32_t minQueryLength) {
        this->minQueryLength = minQueryLength;
    }

    void setQueryLengthStats(int32_t min, double mean, int32_t max);
    
    void setRefs(vector<RefSeq> refs) {
        this->refs = refs;
    }
    
    
    bool addJunction(JunctionPtr j);
    
    /**
     * Appends a new copy of all the junctions in the other junction system to this
     * junction system, without the Bam alignments associated with them.
     * @param other The other junctions system containing junctions to be added to this.
     */
    void append(JunctionSystem& other);
    
    /**
     * Adds any new junctions found from the given alignment to the set managed 
     * by this class
     * @param al The alignment to search for junctions
     * @return Whether a junction was found in this alignment or not
     */
    bool addJunctions(const BamAlignment& al, StrandSpecific strandSpecific) {
        return addJunctions(al, 0, al.getPosition(), strandSpecific);
    }
    
    bool addJunctions(const BamAlignment& al, const size_t startOp, const int32_t offset, StrandSpecific strandSpecific);
    
    
    void findFlankingAlignments(const path& alignmentsFile, StrandSpecific strandSpecific) {
        findFlankingAlignments(alignmentsFile, strandSpecific, false);
    }
    
    void findFlankingAlignments(const path& alignmentsFile, StrandSpecific strandSpecific, bool verbose);
    
    void calcCoverage(const path& alignmentsFile, StrandSpecific strandSpecific);
    
    void calcMultipleMappingStats(SplicedAlignmentMap& map) {
        calcMultipleMappingStats(map, false);
    }
    
    void calcMultipleMappingStats(SplicedAlignmentMap& map, bool verbose);
    
    void calcJunctionStats() {
        calcJunctionStats(false);
    }
    
    void calcJunctionStats(bool verbose);
    
    void calcTrimmedOverhangScore();
    
    void sort();
    
    void saveAll(const path& outputPrefix);
    
    void outputDescription(std::ostream &strm);
    
    friend std::ostream& operator<<(std::ostream &strm, const JunctionSystem& js) {
        
        strm << "index\t" << Junction::junctionOutputHeader() << endl;
        
        uint64_t i = 0;
        for(const auto& j : js.junctionList) {
            strm << i++ << "\t" << *j << endl;
        }
        
        return strm;
    }
    
    void outputJunctionGFF(std::ostream &strm);
    
    void outputIntronGFF(std::ostream &strm);
    
    void outputGTF(std::ostream &strm) {
        
    }
    
    void outputBED(string& path, CanonicalSS type);
    
    void outputBED(std::ostream &strm, CanonicalSS type);
    
    void load(const path& junctionTabFile);
    
    JunctionPtr getJunctionAt(uint32_t index) const {
        return this->junctionList[index];
    }
    
    JunctionPtr getJunction(Intron& intron) const;
    
};
}
