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
#include <vector>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <api/BamAlignment.h>

#include "location.hpp"

using std::endl;
using std::min;
using std::max;
using std::string;
using std::size_t;
using std::vector;

using boost::lexical_cast;
using boost::shared_ptr;

using namespace BamTools;

using portculis::Location;

namespace portculis {    

const uint16_t MAP_QUALITY_THRESHOLD = 30;    
    
class Junction {
private:
    
    // **** Properties that describe where the junction is ****
    shared_ptr<Location> intron;
    vector<BamAlignment> junctionAlignments;
    
    
    // **** Junction metrics ****
                                                // Metric 1 (nbReads) derived from size of junction alignment vector
    bool donorAndAcceptorMotif;                 // Metric 2 calculated from 
                                                // Metric 3 (intron size) calculated via location properties
    int32_t maxMinAnchor;                       // Metric 4
    int32_t diffAnchor;                         // Metric 5
    double entropy;                             // Metric 6
    uint32_t nbDistinctAnchors;                 // Metric 7
    uint32_t nbDistinctAlignments;              // Metric 8
    uint32_t nbReliableAlignments;                // Metric 9
    uint32_t nbUpstreamFlankingAlignments;      // Metric 10
    uint32_t nbDownstreamFlankingAlignments;    // Metric 11
    
    
    
    // **** Additional properties ****
    
    int32_t leftFlankStart;
    int32_t rightFlankEnd;
    string da1, da2;                    // These store the nucleotides found at the predicted donor / acceptor sites in the intron
    
    
    
    void init(int32_t _leftFlankStart, int32_t _rightFlankEnd) {
        
        leftFlankStart = _leftFlankStart;
        rightFlankEnd = _rightFlankEnd;
        donorAndAcceptorMotif = false;
        maxMinAnchor = minAnchor(_leftFlankStart, _rightFlankEnd);
        diffAnchor = 0;
        entropy = 0;
        nbDistinctAnchors = 0;
        nbDistinctAlignments = 0;
        nbReliableAlignments = 0;
        nbUpstreamFlankingAlignments = 0;
        nbDownstreamFlankingAlignments = 0;
    }
    
    int32_t minAnchor(int32_t otherStart, int32_t otherEnd) {
        
        if (intron != NULL) {
            int32_t lAnchor = intron->start - 1 - otherStart;
            int32_t rAnchor = otherEnd - intron->end - 1;

            return min(lAnchor, rAnchor);
        }
        
        return 0;
    }
    
    
protected:
    
    
    
public:
    
    // **** Constructors ****
    
    Junction() {
        init(0, 0);
    }
    
    Junction(shared_ptr<Location> _location, int32_t _leftFlankStart, int32_t _rightFlankEnd) :
        intron(_location) {
        init(_leftFlankStart, _rightFlankEnd);
    }
    
    // **** Destructor ****
    virtual ~Junction() {
        
    }
    
   
    shared_ptr<Location> getLocation() const {
        return intron;
    }

    void setLocation(shared_ptr<Location> location) {
        this->intron = location;
    }

    
    
    void addJunctionAlignment(const BamAlignment& al) {
        this->junctionAlignments.push_back(al);                        
    }
    
    void setDonorAndAcceptorMotif(bool donorAndAcceptorMotif) {
        this->donorAndAcceptorMotif = donorAndAcceptorMotif;
    }
    
    bool setDonorAndAcceptorMotif(string seq1, string seq2) {
        this->da1 = seq1;
        this->da2 = seq2;
        this->donorAndAcceptorMotif = validDonorAcceptor(seq1, seq2);
        return this->donorAndAcceptorMotif;
    }
    
    void setFlankingAlignmentCounts(uint32_t nbUpstream, uint32_t nbDownstream) {
        this->nbUpstreamFlankingAlignments = nbUpstream;
        this->nbDownstreamFlankingAlignments = nbDownstream;
    }

    
    /**
     * Extends the flanking regions of the junction.  Also updates any relevant 
     * metrics.
     * @param otherStart The alternative start position of the left flank
     * @param otherEnd The alternative end position of the right flank
     */
    void extendFlanks(int32_t otherStart, int32_t otherEnd) {        
        
        leftFlankStart = min(leftFlankStart, otherStart);
        rightFlankEnd = max(rightFlankEnd, otherEnd);
        
        int32_t otherMinAnchor = minAnchor(otherStart, otherEnd);
        
        maxMinAnchor = max(maxMinAnchor, otherMinAnchor);
    }
        
    
    /**
     * Call this method to recalculate all junction metrics based on the current location
     * and alignment information present in this junction
     */
    void calcAllMetrics() {
       
        calcAnchorStats();      // Metrics 5 and 7
        calcEntropy();          // Metric 6
        calcAlignmentStats();   // Metrics 8 and 9
    }
    
    /**
     * Metric 5 and 7: Diff Anchor and # Distinct Anchors
     * @return 
     */
    void calcAnchorStats() {
        
        int32_t lEnd = intron->start;
        int32_t rStart = intron->end;
        int32_t minLeftSize = intron->start - leftFlankStart; 
        int32_t minRightSize = rightFlankEnd - intron->end; 
        int32_t maxLeftSize = 0, maxRightSize = 0;
        size_t nbAlignments = junctionAlignments.size();
        uint32_t nbDistinctLeftAnchors = 0, nbDistinctRightAnchors = 0;
        int32_t lastLStart = -1, lastREnd = -1;
                
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            
            int32_t lStart = max(ba.Position, leftFlankStart);
            int32_t rEnd = min(ba.Position + (int32_t)ba.AlignedBases.size(), rightFlankEnd);
            int32_t leftSize = lEnd - lStart;
            int32_t rightSize = rEnd - rStart;
            minLeftSize = min(minLeftSize, leftSize);
            minRightSize = min(minRightSize, rightSize);
            maxLeftSize = max(maxLeftSize, leftSize);
            maxRightSize = max(maxRightSize, rightSize);
            
            if (lStart != lastLStart) {
                nbDistinctLeftAnchors++;
                lastLStart = lStart;
            }
            if (rEnd != lastREnd) {
                nbDistinctRightAnchors++;
                lastREnd = rEnd;
            }
        }
        
        int32_t diffLeftSize = maxLeftSize - minLeftSize;
        int32_t diffRightSize = maxRightSize - minRightSize;
        diffAnchor = min(diffLeftSize, diffRightSize);
        nbDistinctAnchors = min(nbDistinctLeftAnchors, nbDistinctRightAnchors);
    }
    
    
    /**
     * Metric 6: Entropy (definition from "Graveley et al, The developmental 
     * transcriptome of Drosophila melanogaster, Nature, 2011")
     * 
     * We measured the entropy of the reads that mapped to the splice junction. 
     * The entropy score is a function of both the total number of reads 
     * that map to a given junction and the number of different offsets to which 
     * those reads map and the number that map at each offset. Thus, junctions 
     * with multiple reads mapping at each of the possible windows across the junction 
     * will be assigned a higher entropy score, than junctions where many reads 
     * map to only one or two positions. For this analysis we required that a junction 
     * have an entropy score of two or greater in at least two biological samples 
     * for junctions with canonical splice sites, and an entropy score of three 
     * or greater in at least three biological samples for junctions with non-canonical 
     * splice sites. Entropy was calculated using the following equations: 
     * 
     * pi = reads at offset i / total reads to junction window 
     * 
     * Entopy = - sumi(pi * log(pi) / log2) 
     * 
     * @return The entropy of this junction
     */
    double calcEntropy() {
        size_t nbJunctionAlignments = junctionAlignments.size();
        
        if (nbJunctionAlignments <= 1)
            return 0;
        
        double sum = 0.0;
        
        int32_t lastOffset = max(junctionAlignments[0].Position, leftFlankStart);
        uint32_t readsAtOffset = 0;
        
        for(size_t i = 0; i < nbJunctionAlignments; i++) {
            
            const BamAlignment* ba = &(junctionAlignments[i]);
            int32_t pos = max(ba->Position, leftFlankStart);
            
            readsAtOffset++;
            
            if (pos != lastOffset || i == nbJunctionAlignments - 1) {
                double pI = (double)readsAtOffset / (double)nbJunctionAlignments;
                sum += pI * log2(pI);
                lastOffset = pos;
                readsAtOffset = 0;
            }
        }
        
        entropy = -sum;        
        return entropy;
    }
     
    
    /**
     * Metric 8, 9: # Distinct Alignments, # Unique Alignments
     * @return 
     */
    void calcAlignmentStats() {
        
        int32_t lastStart = -1, lastEnd = -1;
        
        nbDistinctAlignments = 0;
        nbReliableAlignments = 0;
        
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            
            int32_t start = ba.Position;
            int32_t end = ba.Position + ba.AlignedBases.size();
            
            if (start != lastStart || end != lastEnd ) {
                nbDistinctAlignments++;                
                lastStart = start;
                lastEnd = end;
            }

            // TODO Suspect this is wrong!!
            // This doesn't seem intuitive but this is how they recommend finding
            // "reliable" (i.e. unique) alignments in samtools.  They do this
            // because apparently "uniqueness" is not a well defined concept.
            if (ba.MapQuality >= MAP_QUALITY_THRESHOLD) {
                nbReliableAlignments++;
            }
        }        
    }
    
    double coverage(int32_t a, int32_t b) {
        
        double multiplier = 1.0 / (b - a);
        uint32_t readCount = 0.0;
        for (int32_t i = a; i <= b; i++) {
            readCount += 0; // Number of reads found at this position
        }
        return multiplier * (double)readCount;
    }
    
    void calcCoverage() {
        
        int32_t readLength = 0;
        uint32_t donorCoverage = coverage(intron->start - 2 * readLength, intron->start - readLength) -
                coverage(intron->start - readLength, intron->start);
        uint32_t acceptorCoverage = coverage(intron->end, intron->end + readLength) -
                coverage(intron->end - readLength, intron->end);
        
        uint32_t coverage = donorCoverage + acceptorCoverage;
    }
    
    
    // **** Core property getters ****
    
    shared_ptr<Location> getIntron() const {
        return intron;
    }

    int32_t getLeftFlankStart() const {
        return leftFlankStart;
    }

    int32_t getRightFlankEnd() const {
        return rightFlankEnd;
    }

    
    
    // **** Metric getters ****
    
    /**
     * Metric 1: The number of alignments directly supporting this junction
     * @return 
     */
    size_t getNbJunctionAlignments() const {
        return this->junctionAlignments.size();
    }
    
    /**
     * Metric 2: Whether or not there is a donor and acceptor motif at the two base
     * pairs at the start and end of the junction / intron
     * @return 
     */
    bool hasDonorAndAcceptorMotif() const {
        return this->donorAndAcceptorMotif;
    }
    
    /**
     * Metric 3: The intron size
     * @return 
     */
    int32_t getIntronSize() const {
        return intron != NULL ? intron->size() : 0;
    }
    
    /**
     * Metric 4: The maximum anchor distance from the shortest side of each supporting 
     * alignment
     * @return
     */
    int32_t getMaxMinAnchor() const {
        return this->maxMinAnchor;
    }
    
    /**
     * Metric 5: Diff Anchor
     * @return 
     */
    int32_t getDiffAnchor() const {        
        return this->diffAnchor;
    }
    
    /**
     * Metric 6: Entropy (definition from "Graveley et al, The developmental 
     * transcriptome of Drosophila melanogaster, Nature, 2011")
     * 
     * We measured the entropy of the reads that mapped to the splice junction. 
     * The entropy score is a function of both the total number of reads 
     * that map to a given junction and the number of different offsets to which 
     * those reads map and the number that map at each offset. Thus, junctions 
     * with multiple reads mapping at each of the possible windows across the junction 
     * will be assigned a higher entropy score, than junctions where many reads 
     * map to only one or two positions. For this analysis we required that a junction 
     * have an entropy score of two or greater in at least two biological samples 
     * for junctions with canonical splice sites, and an entropy score of three 
     * or greater in at least three biological samples for junctions with non-canonical 
     * splice sites. Entropy was calculated using the following equations: 
     * 
     * pi = reads at offset i / total reads to junction window 
     * 
     * Entopy = - sumi(pi * log(pi) / log2) 
     * 
     * @return 
     */
    double getEntropy() const {        
        return this->entropy;
    }
    
    
    /**
     * Metric 7
     * @return 
     */
    uint32_t getNbDistinctAnchors() const {
        return nbDistinctAnchors;
    }
    
    /**
     * Metric 8
     * @return 
     */
    uint32_t getNbDistinctAlignments() const {
        return nbDistinctAlignments;
    }

    /**
     * Metric 9
     * @return 
     */
    uint32_t getNbReliableAlignments() const {
        return nbReliableAlignments;
    }

    /**
     * Metric 10: The number of upstream non-spliced supporting reads
     * @return 
     */
    uint32_t getNbUpstreamFlankingAlignments() const {
        return nbUpstreamFlankingAlignments;
    }
    
    /**
     * Metric 11: The number of downstream non-spliced supporting reads
     * @return 
     */
    uint32_t getNbDownstreamFlankingAlignments() const {
        return nbDownstreamFlankingAlignments;
    }

    
    
    // **** Output methods ****
    
    /**
     * Complete human readable description of this junction
     * @param strm
     */
    void outputDescription(std::ostream &strm) {
        strm << "Location: ";
        
        if (intron != NULL) {
            intron->outputDescription(strm);
        }
        else {
            strm << "No location set";
        }
        
        strm << endl
             << "Flank limits: (" << leftFlankStart << ", " << rightFlankEnd << ")" << endl
             << "Junction Metrics:" << endl
             << "1:  # Junction Alignments: " << getNbJunctionAlignments() << endl
             << "2:  Has Donor + Acceptor Motif: " << donorAndAcceptorMotif << "; Sequences: (" << da1 << " " << da2 << ")" << endl
             << "3:  Intron Size: " << getIntronSize() << endl
             << "4:  MaxMinAnchor: " << maxMinAnchor << endl
             << "5:  DiffAnchor: " << diffAnchor << endl
             << "6:  Entropy: " << entropy << endl
             << "7:  # Distinct Anchors: " << nbDistinctAnchors << endl
             << "8:  # Distinct Alignments: " << nbDistinctAlignments << endl
             << "9:  # Reliable (MP >=" << MAP_QUALITY_THRESHOLD << ") Alignments: " << nbReliableAlignments << endl
             << "10: # Upstream Non-Spliced Alignments: " << nbUpstreamFlankingAlignments << endl
             << "11: # Downstream Non-Spliced Alignments: " << nbDownstreamFlankingAlignments << endl
             << endl;
                
    }
    
    
    /**
     * Represents this junction as a table row
     * @param seq1
     * @param seq2
     * @return 
     */
    friend std::ostream& operator<<(std::ostream &strm, const Junction& j) {
        return strm << *(j.intron) << "\t"
                    << j.leftFlankStart << "\t"
                    << j.rightFlankEnd << "\t"
                    << j.getNbJunctionAlignments() << "\t"
                    << j.donorAndAcceptorMotif << "\t"
                    << j.getIntronSize() << "\t"
                    << j.maxMinAnchor << "\t"
                    << j.diffAnchor << "\t"
                    << j.entropy << "\t"
                    << j.nbDistinctAnchors << "\t"
                    << j.nbDistinctAlignments << "\t"
                    << j.nbReliableAlignments << "\t"
                    << j.nbUpstreamFlankingAlignments << "\t"
                    << j.nbDownstreamFlankingAlignments;
    }
    
        
    // **** Static methods ****

    /**
     * Tests whether the two strings could represent valid donor and acceptor sites
     * @param seq1
     * @param seq2
     * @return 
     */
    static bool validDonorAcceptor(string seq1, string seq2) {
        return ((seq1.size() == 2 && seq2.size() == 2) &&
                ((seq1[0] == 'G' && seq1[1] == 'T' && seq2[0] == 'A' && seq2[1] == 'G') ||
                 (seq1[0] == 'C' && seq1[1] == 'T' && seq2[0] == 'A' && seq2[1] == 'C')));
    }
    
    /**
     * Header for table output
     * @return 
     */
    static string junctionOutputHeader() {
        return string(Location::locationOutputHeader()) + string("\tleft\tright\tM1\tM2\tM3\tM4\tM5\tM6\tM7\tM8\tM9\tM10\tM11"); 
    }

};

}
