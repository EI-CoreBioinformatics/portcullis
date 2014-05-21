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

    
class Junction {
private:
    
    // **** Properties that describe where the junction is ****
    shared_ptr<Location> intron;
    vector<BamAlignment> junctionAlignments;
    vector<BamAlignment> leftFlankingAlignments;
    vector<BamAlignment> rightFlankingAlignments;
    
    
    // **** Junction metrics ****
                                        // Metric 1 (nbReads) derived from size of junction alignment vector
    bool donorAndAcceptorMotif;         // Metric 2 calculated from 
                                        // Metric 3 (intron size) calculated via location properties
    int32_t maxMinAnchor;               // Metric 4
    int32_t diffAnchor;                 // Metric 5
    double entropy;                     // Metric 6
    
    
    // **** Additional stuff ****
    
    int32_t leftFlankStart;
    int32_t rightFlankEnd;
    string da1, da2;                    // These store the nucleotides found at the predicted donor / acceptor sites in the intron
    
    
    
    void init(int32_t _leftFlankStart, int32_t _rightFlankEnd) {
        
        donorAndAcceptorMotif = false;
        leftFlankStart = _leftFlankStart;
        rightFlankEnd = _rightFlankEnd;
        maxMinAnchor = minAnchor(_leftFlankStart, _rightFlankEnd);
        diffAnchor = 0;
        entropy = 0;
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
       
        calcDiffAnchor();   // Metric 5
        calcEntropy();      // Metric 6
    }
    
    /**
     * Metric 5: Diff Anchor
     * @return 
     */
    int32_t calcDiffAnchor() {
        
        int32_t minLeft = intron->start - leftFlankStart; 
        int32_t minRight = rightFlankEnd - intron->end; 
        int32_t maxLeft = 0, maxRight = 0;
        
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            
            int32_t end = ba.Position + ba.Length;
            int32_t left = intron->start - ba.Position;
            int32_t right = max(end, rightFlankEnd) - intron->end;
            minLeft = min(minLeft, left);
            minRight = min(minRight, right);
            maxLeft = max(maxLeft, left);
            maxRight = max(maxRight, right);
        }
        
        int32_t diffLeft = maxLeft - minLeft;
        int32_t diffRight = maxRight - minRight;
        diffAnchor = min(diffLeft, diffRight);
        return diffAnchor;
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
        
        int32_t lastOffset = junctionAlignments[0].Position;
        uint32_t readsAtOffset = 0;
        
        for(size_t i = 0; i < nbJunctionAlignments; i++) {
            
            const BamAlignment* ba = &(junctionAlignments[i]);
            int32_t pos = ba->Position;
            
            if (pos == lastOffset) {
                readsAtOffset++;
            }
            
            if (pos != lastOffset || i == nbJunctionAlignments - 1) {
                double pI = (double)readsAtOffset / (double)getNbJunctionAlignments();
                sum += pI * log2(pI);
                lastOffset = pos;
                readsAtOffset = 0;
            }
        }
        
        entropy = -sum;        
        return entropy;
    }
    
    
    
    // **** Metric getters ****
    
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
        return intron != NULL ? intron->size() : 0;
    }
    
    /**
     * Metric 4: The maximum anchor distance from the shortest side of each supporting 
     * alignment
     * @return
     */
    int32_t getMaxMinAnchor() const {
        return maxMinAnchor;
    }
    
    /**
     * Metric 5: Diff Anchor
     * @return 
     */
    int32_t getDiffAnchor() const {
        
        int32_t minLeft = intron->start - leftFlankStart; 
        int32_t minRight = rightFlankEnd - intron->end; 
        int32_t maxLeft = 0, maxRight = 0;
        
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            
            int32_t end = ba.Position + ba.Length;
            int32_t left = intron->start - ba.Position;
            int32_t right = max(end, rightFlankEnd) - intron->end;
            minLeft = min(minLeft, left);
            minRight = min(minRight, right);
            maxLeft = max(maxLeft, left);
            maxRight = max(maxRight, right);
        }
        
        int32_t diffLeft = maxLeft - minLeft;
        int32_t diffRight = maxRight - minRight;
        return min(diffLeft, diffRight);
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
        return entropy;
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
             << "10: # Upstream Non-Spliced Alignments: " << getNbUpstream() << endl
             << "11: # Downstream Non-Spliced Alignments: " << getNbDownstream() << endl
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
                    << j.getNbUpstream() << "\t"
                    << j.getNbDownstream();
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
        return string(Location::locationOutputHeader()) + string("left\tright\tM1\tM2\tM3\tM4\tM5\tM6\tM10\tM11"); 
    }

};

}
