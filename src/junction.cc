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

#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
using std::endl;
using std::min;
using std::max;
using std::string;
using std::size_t;
using std::vector;
using std::shared_ptr;

typedef std::unordered_map<size_t, uint16_t> SplicedAlignmentMap;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

#include "intron.hpp"
#include "samtools_helper.hpp"
#include "seq_utils.hpp"
#include "prepare.hpp"
using portcullis::GenomeMapper;
using portcullis::Intron;
using portcullis::Strand;
using portcullis::SamtoolsHelper;
using portcullis::CigarOp;

#include "junction.hpp"

/**
 * Tests whether the two strings could represent valid donor and acceptor sites
 * for this junction
 * @param seq1
 * @param seq2
 * @return 
 */
portcullis::CanonicalSS portcullis::Junction::hasCanonicalSpliceSites(const string& seq1, const string& seq2) {

    if (intron == nullptr || seq1.size() != 2 || seq2.size() != 2)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't test for valid donor / acceptor when either string are not of length two, or the intron location is not defined")));

    const string seq = string(seq1 + seq2);

    if (seq == CANONICAL_SEQ || seq == CANONICAL_SEQ_RC) {
        return CANONICAL;
    }
    else if (seq == SEMI_CANONICAL_SEQ1 || seq == SEMI_CANONICAL_SEQ1_RC ||
             seq == SEMI_CANONICAL_SEQ2 || seq == SEMI_CANONICAL_SEQ2_RC) {
        return SEMI_CANONICAL;
    }
    else {
        return NO;
    }
}
    
Strand portcullis::Junction::predictedStrandFromSpliceSites(const string& seq1, const string& seq2) {

    if (seq1.size() != 2 || seq2.size() != 2)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't test donor / acceptor when either string are not of length two")));

    const string seq = string(seq1 + seq2);

    if (seq == CANONICAL_SEQ) {
        return POSITIVE;
    }
    else if (seq == CANONICAL_SEQ_RC) {
        return NEGATIVE;
    }
    // Check for these after canonicals for performance reasons
    else if (seq == SEMI_CANONICAL_SEQ1 || seq == SEMI_CANONICAL_SEQ2) {
        return POSITIVE;
    }
    else if (seq == SEMI_CANONICAL_SEQ1_RC || seq == SEMI_CANONICAL_SEQ2_RC) {
        return NEGATIVE;
    }
    else {
        return UNKNOWN;
    }
}
    

portcullis::Junction::Junction(shared_ptr<Intron> _location, int32_t _leftFlankStart, int32_t _rightFlankEnd) :
    intron(_location) {
    leftFlankStart = _leftFlankStart;
    rightFlankEnd = _rightFlankEnd;
    canonicalSpliceSites = NO;
    maxMinAnchor = intron->minAnchorLength(_leftFlankStart, _rightFlankEnd);
    diffAnchor = 0;
    entropy = 0;
    nbDistinctAnchors = 0;
    nbJunctionAlignments = 0;
    nbDistinctAlignments = 0;
    nbReliableAlignments = 0;
    nbUpstreamFlankingAlignments = 0;
    nbDownstreamFlankingAlignments = 0;
    leftAncSize = 0;
    rightAncSize = 0;
    maxMMES = 0;
    hammingDistance5p = -1;
    hammingDistance3p = -1;
    coverage = 0.0;
    uniqueJunction = false;
    primaryJunction = false;
    multipleMappingScore = 0.0;
    nbMismatches = 0;
    nbMultipleSplicedReads = 0;

    predictedStrand = UNKNOWN;
    nbUpstreamJunctions = 0;
    nbDownstreamJunctions = 0;
}
   
/**
 * Copy constructor, with option to copy bam alignments associated with the junction
 * @param j The other junction to deep copy into this
 * @param withAlignments Whether to copy over the alignments or not
 */
portcullis::Junction::Junction(const Junction& j, bool withAlignments) {
    intron = make_shared<Intron>(*(j.intron));        
    leftFlankStart = j.leftFlankStart;
    rightFlankEnd = j.rightFlankEnd;
    canonicalSpliceSites = j.canonicalSpliceSites;
    maxMinAnchor = j.maxMinAnchor;
    diffAnchor = j.diffAnchor;
    entropy = j.entropy;
    nbDistinctAnchors = j.nbDistinctAnchors;
    nbJunctionAlignments = j.nbJunctionAlignments;
    nbDistinctAlignments = j.nbDistinctAlignments;
    nbReliableAlignments = j.nbReliableAlignments;
    leftAncSize = j.leftAncSize;
    rightAncSize = j.rightAncSize;
    nbUpstreamFlankingAlignments = j.nbUpstreamFlankingAlignments;
    nbDownstreamFlankingAlignments = j.nbDownstreamFlankingAlignments;
    maxMMES = j.maxMMES;
    hammingDistance5p = j.hammingDistance5p;
    hammingDistance3p = j.hammingDistance3p;
    coverage = j.coverage;
    uniqueJunction = j.uniqueJunction;
    primaryJunction = j.primaryJunction;
    multipleMappingScore = j.multipleMappingScore;
    nbMismatches = j.nbMismatches;
    nbMultipleSplicedReads = j.nbMultipleSplicedReads;

    predictedStrand = j.predictedStrand;
    nbUpstreamJunctions = j.nbUpstreamJunctions;
    nbDownstreamJunctions = j.nbDownstreamJunctions;

    if (withAlignments) {

        for(shared_ptr<BamAlignment> ba : j.junctionAlignments) {
            this->junctionAlignments.push_back(ba);
        }

        for(size_t baNameCode : j.junctionAlignmentNames) {
            this->junctionAlignmentNames.push_back(baNameCode);
        }
    }
}
    
// **** Destructor ****
portcullis::Junction::~Junction() {
    junctionAlignments.clear();
    junctionAlignmentNames.clear();
}
    
void portcullis::Junction::clearAlignments() {
    junctionAlignments.clear();  
}
       
    
void portcullis::Junction::addJunctionAlignment(BamAlignmentPtr al) {

    // Make sure we take a proper copy of this alignment for safe storage
    this->junctionAlignments.push_back(make_shared<BamAlignment>(*al));

    // Calculate a hash of the alignment name
    size_t code = std::hash<std::string>()(al->deriveName());

    this->junctionAlignmentNames.push_back(code);
    this->nbJunctionAlignments = this->junctionAlignments.size();

    if (al->getNbJunctionsInRead() > 1) {
        this->nbMultipleSplicedReads++;
    }
}
    
    
portcullis::CanonicalSS portcullis::Junction::setDonorAndAcceptorMotif(string seq1, string seq2) {
    this->da1 = seq1;
    this->da2 = seq2;
    this->canonicalSpliceSites = hasCanonicalSpliceSites(seq1, seq2);
    this->predictedStrand = predictedStrandFromSpliceSites(seq1, seq2);
    return this->canonicalSpliceSites;
}
    
void portcullis::Junction::setFlankingAlignmentCounts(uint32_t nbUpstream, uint32_t nbDownstream) {
    this->nbUpstreamFlankingAlignments = nbUpstream;
    this->nbDownstreamFlankingAlignments = nbDownstream;
}
    

    
/**
 * Extends the flanking regions of the junction.  Also updates any relevant 
 * metrics.
 * @param otherStart The alternative start position of the left flank
 * @param otherEnd The alternative end position of the right flank
 */
void portcullis::Junction::extendFlanks(int32_t otherStart, int32_t otherEnd) {        

    leftFlankStart = min(leftFlankStart, otherStart);
    rightFlankEnd = max(rightFlankEnd, otherEnd);

    int32_t otherMinAnchor = intron->minAnchorLength(otherStart, otherEnd);

    maxMinAnchor = max(maxMinAnchor, otherMinAnchor);
}
        
    
portcullis::CanonicalSS portcullis::Junction::processJunctionWindow(GenomeMapper* genomeMapper) {

    if (intron == nullptr) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find genomic sequence for this junction as no intron is defined")));

    int32_t refid = intron->ref.index;

    // Just access the whole junction region
    int seqLen = -1;
    string region(genomeMapper->fetchBases(intron->ref.name.c_str(), leftFlankStart, rightFlankEnd, &seqLen));
    boost::to_upper(region);    // Removes any lowercase bases representing repeats
    if (seqLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find genomic region for junction")));

    if (seqLen != rightFlankEnd - leftFlankStart + 1)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence is not of the expected length.") +
                "\nSequence Name: " + intron->ref.name + 
                "\nSequence Length: " + lexical_cast<string>(seqLen) + 
                "\nLeft flank start: " + lexical_cast<string>(leftFlankStart) + 
                "\nRight flank end: " + lexical_cast<string>(rightFlankEnd) + "\n"));


    // Process the predicted donor / acceptor regions and update junction
    string daSeq1 = region.substr(intron->start - leftFlankStart, 2);
    string daSeq2 = region.substr(intron->end - 2 - leftFlankStart, 2);
    CanonicalSS validDA = setDonorAndAcceptorMotif(daSeq1, daSeq2);

    // Create strings for hamming analysis
    calcHammingScores(region);

    return validDA;
}
    
void portcullis::Junction::processJunctionVicinity(BamReader& reader, int32_t refLength, int32_t meanQueryLength, int32_t maxQueryLength, StrandSpecific strandSpecific) {

    int32_t refId = intron->ref.index;
    Strand strand = intron->strand;

    uint32_t nbLeftFlankingAlignments = 0, nbRightFlankingAlignments = 0;

    int32_t regionStart = leftFlankStart - maxQueryLength - 1; regionStart < 0 ? 0 : regionStart;
    int32_t regionEnd = rightFlankEnd + maxQueryLength + 1; regionEnd >= refLength ? refLength - 1 : regionEnd;

    // Focus only on the (expanded... to be safe...) region of interest
    reader.setRegion(refId, regionStart, regionEnd);

    while(reader.next()) {

        BamAlignmentPtr ba = reader.current();
        int32_t pos = ba->getPosition();

        //TODO: Should we consider strand specific reads differently here?

        // Look for left flanking alignments
        if (    intron->start > pos && 
                leftFlankStart <= pos + ba->calcNbAlignedBases()) {
            nbLeftFlankingAlignments++;
        }

        // Look for right flanking alignments
        if (    rightFlankEnd >= pos && 
                intron->end < pos) {
            nbRightFlankingAlignments++;
        }
    }

    this->setFlankingAlignmentCounts(nbLeftFlankingAlignments, nbRightFlankingAlignments);
}
    
    
void portcullis::Junction::calcMetrics() {

    calcAnchorStats();      // Metrics 5 and 7
    calcEntropy();          // Metric 6
    calcAlignmentStats();   // Metrics 8, 9 and 19
    calcMaxMMES();          // Metric 12
}
    
    
/**
 * Metric 5 and 7: Diff Anchor and # Distinct Anchors
 * @return 
 */
void portcullis::Junction::calcAnchorStats() {

    int32_t lEnd = intron->start;
    int32_t rStart = intron->end;
    int32_t minLeftSize = intron->start - leftFlankStart; 
    int32_t minRightSize = rightFlankEnd - intron->end; 
    int32_t maxLeftSize = 0, maxRightSize = 0;
    size_t nbAlignments = junctionAlignments.size();
    uint32_t nbDistinctLeftAnchors = 0, nbDistinctRightAnchors = 0;
    int32_t lastLStart = -1, lastREnd = -1;

    for(BamAlignmentPtr ba : junctionAlignments) {

        const int32_t pos = ba->getPosition();
        const int32_t alignedBases = (int32_t)ba->calcNbAlignedBases();

        int32_t leftSize = ba->alignedBasesBetween(pos, lEnd);
        int32_t rightSize = ba->alignedBasesBetween(rStart, pos + alignedBases);

        maxLeftSize = max(maxLeftSize, leftSize);
        minLeftSize = min(minLeftSize, leftSize);
        maxRightSize = max(maxRightSize, rightSize);
        minRightSize = min(minRightSize, rightSize);

        int32_t lStart = pos;
        int32_t rEnd = pos + alignedBases;

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
    leftAncSize = maxLeftSize;
    rightAncSize = maxRightSize;            
    diffAnchor = min(diffLeftSize, diffRightSize);
    nbDistinctAnchors = min(nbDistinctLeftAnchors, nbDistinctRightAnchors);
}
    
    
/**
 * Metric 6: Entropy (definition from "Graveley et al, The developmental 
 * transcriptome of Drosophila melanogaster, Nature, 2011")
 * 
 * Calculates the entropy score for this junction.  Higher entropy is generally
 * more indicative of a genuine junction than a lower score.
 *
 * This version of the function collects BamAlignment start positions and passes
 * them to an overloaded version of this function.
 * 
 * @return The entropy of this junction
 */
double portcullis::Junction::calcEntropy() {

    vector<int32_t> junctionPositions;

    for(BamAlignmentPtr ba : junctionAlignments) {
        junctionPositions.push_back(ba->getPosition());
    }

    return calcEntropy(junctionPositions);        
}
    
/**
 * Metric 6: Entropy (definition from "Graveley et al, The developmental 
 * transcriptome of Drosophila melanogaster, Nature, 2011")
 *
 * Calculates the entropy score for this junction.  Higher entropy is generally
 * more indicative of a genuine junction than a lower score.
 * 
 * This overridden version of the calcEntropy method, allows us to calculate
 * the entropy, without using BamAlignment objects that model the junction.
 * All we need are the alignment start positions instead.
 * 
 * We measured the entropy of the reads that mapped to the splice junction. 
 * The entropy score is a function of both the total number of reads 
 * that map to a given junction and the number of different offsets to which 
 * those reads map and the number that map at each offset. Thus, junctions 
 * with multiple reads mapping at each of the possible windows across the junction 
 * will be assigned a higher entropy score, than junctions where many reads 
 * map to only one or two positions. 
 * 
 * Entropy was calculated using the following equations: 
 * 
 * p_i = nb_reads_at_offset_i / total_reads_in_junction_window 
 * 
 * Entropy = - sum_i(p_i * log(pi) / log2) 
 * 
 * @param junctionPositions start index positions of each junction alignment
 * 
 * @return The entropy of this junction
 */
double portcullis::Junction::calcEntropy(vector<int32_t>& junctionPositions) {

    size_t nbJunctionAlignments = junctionPositions.size();

    if (nbJunctionAlignments <= 1)
        return 0;

    double sum = 0.0;

    int32_t lastOffset = max(junctionPositions[0], leftFlankStart);
    uint32_t readsAtOffset = 0;

    for(size_t i = 0; i < nbJunctionAlignments; i++) {

        int32_t pos = max(junctionPositions[i], leftFlankStart);

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
 * Metrics: # Distinct Alignments, # Unique/Reliable Alignments, #mismatches
 * @return 
 */
void portcullis::Junction::calcAlignmentStats() {

    int32_t lastStart = -1, lastEnd = -1;

    nbDistinctAlignments = 0;
    nbReliableAlignments = 0;
    nbMismatches = 0;
    nbUpstreamJunctions = 0;
    nbDownstreamJunctions = 0;

    //cout << junctionAlignments.size() << endl;

    for(BamAlignmentPtr ba : junctionAlignments) {

        const int32_t bapos = ba->getPosition();
        const int32_t alignedBases = (int32_t)ba->calcNbAlignedBases();

        int32_t start = bapos;
        int32_t end = bapos + alignedBases;

        if (start != lastStart || end != lastEnd ) {
            nbDistinctAlignments++;                
            lastStart = start;
            lastEnd = end;
        }

        // TODO Suspect this is wrong!!
        // This doesn't seem intuitive but this is how they recommend finding
        // "reliable" (i.e. unique) alignments in samtools.  They do this
        // because apparently "uniqueness" is not a well defined concept.
        if (ba->getMapQuality() >= MAP_QUALITY_THRESHOLD) {
            nbReliableAlignments++;
        }

        uint16_t leftSize = 0;
        uint16_t rightSize = 0;
        uint16_t upjuncs = 0;
        uint16_t downjuncs = 0;
        int32_t pos = start;

        for (CigarOp op : ba->getCigar()) {

            if (CigarOp::opConsumesReference(op.type)) {
                pos += op.length;                    
            }

            if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
                if (pos < intron->start) {
                    upjuncs++;
                }
                else if (pos > intron->end) {
                    downjuncs++;
                }  
            }

            // Looks for a mismatch at any point along the spliced read
            // (even if it's the otherside of another junction in the case of an MSR)
            if (op.type == BAM_CIGAR_DIFF_CHAR) {
                nbMismatches++;
            }
        }

        nbUpstreamJunctions = max(nbUpstreamJunctions, upjuncs);
        nbDownstreamJunctions = max(nbDownstreamJunctions, downjuncs);
    }        
}
    
/**
 * Metric 13 and 14: Calculates the 5' and 3' hamming distances from a genomic
 * region represented by this junction
 * @param junctionSeq The DNA sequence representing this junction on the genome
 */
void portcullis::Junction::calcHammingScores(string& junctionSeq) {

    // Default region length is 10, however this may get modified if we have
    // particularly short anchors
    uint16_t REGION_LENGTH = 10;

    string intron5p;
    string intron3p;
    string anchor5p;
    string anchor3p;

    int32_t diffAnchor5p = intron->start - leftFlankStart;
    int32_t diffAnchor3p = rightFlankEnd - intron->end;

    // Shorten region to look for based on size of the anchors (flanking region).  i.e. if the 
    // smallest anchor is less than the default region length, we shrink the region length
    // to the size of the smallest anchor
    REGION_LENGTH = diffAnchor5p < REGION_LENGTH ? diffAnchor5p : REGION_LENGTH;
    REGION_LENGTH = diffAnchor3p < REGION_LENGTH ? diffAnchor3p : REGION_LENGTH;

    //cout << "region length: " << REGION_LENGTH << endl;

    int32_t intronStartOffset = intron->start - leftFlankStart;
    int32_t anchor3pStartOffset = intron->end - leftFlankStart + 1;

    // The first parts of this if statement are for debugging purposes.  We
    // don't actually expect these to be set for real
    if (intron->size() < REGION_LENGTH) {
        hammingDistance5p = -1;
        hammingDistance3p = -1;
    }
    else if (intronStartOffset - REGION_LENGTH < 0) {
        hammingDistance5p = -2;
        hammingDistance3p = -2;
    }
    else if (anchor3pStartOffset + REGION_LENGTH > size()) {
        hammingDistance5p = -3;
        hammingDistance3p = -3;
    }
    else if (anchor3pStartOffset + REGION_LENGTH > junctionSeq.size()) {
        hammingDistance5p = -4;
        hammingDistance3p = -4;
    }
    else {

        Strand s = intron->strand != UNKNOWN ? intron->strand : predictedStrand;                

        switch(s) {
            case POSITIVE:
                intron5p = junctionSeq.substr(intronStartOffset, REGION_LENGTH);
                intron3p = junctionSeq.substr(anchor3pStartOffset - REGION_LENGTH, REGION_LENGTH);
                anchor5p = junctionSeq.substr(intronStartOffset - REGION_LENGTH, REGION_LENGTH);
                anchor3p = junctionSeq.substr(anchor3pStartOffset, REGION_LENGTH);
                break;
            case NEGATIVE:
                intron5p = SeqUtils::reverseComplement(junctionSeq.substr(anchor3pStartOffset - REGION_LENGTH, REGION_LENGTH));
                intron3p = SeqUtils::reverseComplement(junctionSeq.substr(intronStartOffset, REGION_LENGTH));
                anchor5p = SeqUtils::reverseComplement(junctionSeq.substr(anchor3pStartOffset, REGION_LENGTH));
                anchor3p = SeqUtils::reverseComplement(junctionSeq.substr(intronStartOffset - REGION_LENGTH, REGION_LENGTH));
                break;
            default:
                hammingDistance5p = -5;
                hammingDistance3p = -5;
                return;
        }

        /*cout << "anchor 5': " << anchor5p << endl;
        cout << "intron 5': " << intron5p << endl;
        cout << "intron 3': " << intron3p << endl;
        cout << "anchor 3': " << anchor3p << endl;*/

        hammingDistance5p = SeqUtils::hammingDistance(anchor5p, intron3p);
        hammingDistance3p = SeqUtils::hammingDistance(anchor3p, intron5p);  
    }
}
    
/**
 * Calculates metric 12.  MaxMMES.
 */
void portcullis::Junction::calcMaxMMES() {

    uint16_t maxmmes = 0;

    for(BamAlignmentPtr ba : junctionAlignments) {

        uint16_t leftMM = ba->calcMinimalMatchInCigarDataSubset(leftFlankStart, intron->start);
        uint16_t rightMM = ba->calcMinimalMatchInCigarDataSubset(intron->end, rightFlankEnd);

        uint16_t mmes = min(leftMM, rightMM);

        maxmmes = max(maxmmes, mmes);
    }

    this->maxMMES = maxmmes;
}
    
/**
 * Calculates metric 18.  Multiple mapping score
 */
void portcullis::Junction::calcMultipleMappingScore(SplicedAlignmentMap& map) {

    size_t N = this->getNbJunctionAlignments();

    uint32_t M = 0;
    for(size_t baNameCode : junctionAlignmentNames) {
        M += map[baNameCode];  // Number of multiple splitting patterns
    }

    this->multipleMappingScore = (double)N / (double)M;
}
    
    
double portcullis::Junction::calcCoverage(int32_t a, int32_t b, const vector<uint32_t>& coverageLevels) {

    double multiplier = 1.0 / (b - a);
    uint32_t readCount = 0;

    for (int32_t i = a; i <= b; i++) {

        // Don't do anything stupid!
        if (i >= 0 && i < coverageLevels.size()) {
            readCount += coverageLevels[i];
        }            
    }
    return multiplier * (double)readCount;
}
    
double portcullis::Junction::calcCoverage(const vector<uint32_t>& coverageLevels) {

    const uint32_t REGION_LENGTH = 10;

    int32_t donorStart = intron->start - 2 * REGION_LENGTH;
    int32_t donorMid = intron->start - REGION_LENGTH;
    int32_t donorEnd = intron->start;

    int32_t acceptorStart = intron->end;
    int32_t acceptorMid = intron->end + REGION_LENGTH;
    int32_t acceptorEnd = intron->end + 2 * REGION_LENGTH;

    double donorCoverage = 
            calcCoverage(donorStart, donorMid-1, coverageLevels) -
            calcCoverage(donorMid, donorEnd, coverageLevels);

    double acceptorCoverage = 
            calcCoverage(acceptorMid, acceptorEnd, coverageLevels) -
            calcCoverage(acceptorStart, acceptorMid-1, coverageLevels);

    coverage = donorCoverage + acceptorCoverage; 

    return coverage;
}
    
    
    


    
    // **** Output methods ****
    
/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::outputDescription(std::ostream &strm, string delimiter) {

    if (intron != nullptr) {
        intron->outputDescription(strm, delimiter);
    }
    else {
        strm << "No location set";
    }

    strm << delimiter
         << "Flank limits: (" << leftFlankStart << ", " << rightFlankEnd << ")" << delimiter
         << "Junction Predictions:" << delimiter
         << "P1:  Strand: " << strandToString(predictedStrand) << delimiter
         << "Junction Metrics:" << delimiter
         << "M1:  Canonical?: " << boolalpha << cssToString(canonicalSpliceSites) << "; Sequences: (" << da1 << " " << da2 << ")" << delimiter
         << "M2:  # Junction Alignments: " << getNbJunctionAlignments() << delimiter
         << "M3:  # Distinct Alignments: " << nbDistinctAlignments << delimiter
         << "M4:  # Reliable (MP >=" << MAP_QUALITY_THRESHOLD << ") Alignments: " << nbReliableAlignments << delimiter
         << "M5:  Intron Size: " << getIntronSize() << delimiter
         << "M6:  Left Anchor Size: " << leftAncSize << delimiter
         << "M7:  Right Anchor Size: " << rightAncSize << delimiter
         << "M8:  MaxMinAnchor: " << maxMinAnchor << delimiter
         << "M9:  DiffAnchor: " << diffAnchor << delimiter
         << "M10: # Distinct Anchors: " << nbDistinctAnchors << delimiter
         << "M11: Entropy: " << entropy << delimiter
         << "M12: MaxMMES: " << maxMMES << delimiter
         << "M13: Hamming Distance 5': " << hammingDistance5p << delimiter
         << "M14: Hamming Distance 3': " << hammingDistance3p << delimiter
         << "M15: Coverage: " << coverage << delimiter
         << "M16: Unique Junction: " << boolalpha << uniqueJunction << delimiter
         << "M17: Primary Junction: " << boolalpha << primaryJunction << delimiter
         << "M18: Multiple mapping score: " << multipleMappingScore << delimiter
         << "M19: # mismatches: " << nbMismatches << delimiter
         << "M20: # Multiple Spliced Reads: " << nbMultipleSplicedReads << delimiter
         << "M21: # Upstream Junctions: " << nbUpstreamJunctions << delimiter
         << "M22: # Downstream Junctions: " << nbDownstreamJunctions << delimiter
         << "M23: # Upstream Non-Spliced Alignments: " << nbUpstreamFlankingAlignments << delimiter
         << "M24: # Downstream Non-Spliced Alignments: " << nbDownstreamFlankingAlignments;

}
    
    
/**
 * Complete human readable description of this intron (for augustus hints)
 * @param strm
 */
void portcullis::Junction::outputIntronGFF(std::ostream &strm, uint32_t id) {

    // Use intron strand if known, otherwise use the predicted strand,
    // if predicted strand is also unknown then use "." to indicated unstranded
    const char strand = (intron->strand == UNKNOWN ?
                            predictedStrand == UNKNOWN ?
                                '.' :
                                strandToChar(predictedStrand) :
                            strandToChar(intron->strand));

    string juncId = string("junc_") + lexical_cast<string>(id);

    // Output junction parent
    strm << intron->ref.name << "\t"
         << "portcullis" << "\t"    // source
         << "intron" << "\t"        // type (may change later)
         << intron->start << "\t"   // start
         << intron->end << "\t"     // end
         << "0.0" << "\t"           // No score for the moment
         << strand << "\t"          // strand
         << "." << "\t"             // Just put "." for the phase
         << "mult=" << nbJunctionAlignments << ";"  // Number of times it was seen
         << "src=E";                // Source for augustus
    strm << endl;

}

/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::outputJunctionGFF(std::ostream &strm, uint32_t id) {

    // Use intron strand if known, otherwise use the predicted strand,
    // if predicted strand is also unknown then use "." to indicated unstranded
    const char strand = (intron->strand == UNKNOWN ?
                            predictedStrand == UNKNOWN ?
                                '.' :
                                strandToChar(predictedStrand) :
                            strandToChar(intron->strand));

    string juncId = string("junc_") + lexical_cast<string>(id);

    // Output junction parent
    strm << intron->ref.name << "\t"
         << "portcullis" << "\t"    // source
         << "junction" << "\t"      // type (may change later)
         << leftFlankStart << "\t"  // start
         << rightFlankEnd << "\t"   // end
         << "0.0" << "\t"           // No score for the moment
         << strand << "\t"          // strand
         << "." << "\t"             // Just put "." for the phase
         << "ID=" << juncId << ";"  // ID of the intron
         << "mult=" << nbJunctionAlignments << ";"  // Number of times it was seen
         << "src=E;";                // Source for augustus
    outputDescription(strm, ";");
    strm << endl;

    // Output left exonic region
    strm << intron->ref.name << "\t"
         << "portcullis" << "\t"
         << "partial_exon" << "\t"
         << leftFlankStart << "\t"
         << (intron->start - 1) << "\t"
         << "0.0" << "\t"
         << strand << "\t"
         << "." << "\t"
         << "ID=" << juncId << "_left" << ";"
         << "Parent=" << juncId << endl;

    // Output right exonic region
    strm << intron->ref.name << "\t"
         << "portcullis" << "\t"
         << "partial_exon" << "\t"
         << (intron->end + 1) << "\t"
         << rightFlankEnd << "\t"
         << "0.0" << "\t"
         << strand << "\t"
         << "." << "\t"
         << "ID=" << juncId << "_right" << ";"
         << "Parent=" << juncId << endl;
}
    
/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::outputBED(std::ostream &strm, uint32_t id) {

    // Use intron strand if known, otherwise use the predicted strand,
    // if predicted strand is also unknown then use "." to indicated unstranded
    const char strand = (intron->strand == UNKNOWN ?
                            predictedStrand == UNKNOWN ?
                                '.' :
                                strandToChar(predictedStrand) :
                            strandToChar(intron->strand));

    string juncId = string("junc_") + lexical_cast<string>(id);

    int32_t sz1 = intron->start - leftFlankStart;
    int32_t sz2 = rightFlankEnd - intron->end;
    string blockSizes = lexical_cast<string>(sz1) + "," + lexical_cast<string>(sz2);
    string blockStarts = lexical_cast<string>(0) + "," + lexical_cast<string>(intron->end - leftFlankStart);

    // Output junction parent
    strm << intron->ref.name << "\t"         // chrom
         << leftFlankStart << "\t"  // chromstart
         << rightFlankEnd << "\t"   // chromend
         << juncId << "\t"          // name
         << this->getNbJunctionAlignments() << "\t"           // Use the depth as the score for the moment
         << strand << "\t"          // strand
         << intron->start << "\t"   // thickstart
         << intron->end << "\t"     // thickend
         << "255,0,0" << "\t"       // Just use red for the moment
         << "2" << "\t"             // 2 blocks: Left and right block
         << blockSizes << "\t"
         << blockStarts << endl;        
}
    
    
    
// **** Static methods ****

/**
 * Header for table output
 * @return 
 */
string portcullis::Junction::junctionOutputHeader() {
    return string(Intron::locationOutputHeader()) + "\tleft\tright\tss1\tss2\t" + 
            boost::algorithm::join(PREDICTION_NAMES, "\t") + "\t" +
            boost::algorithm::join(METRIC_NAMES, "\t");
}

shared_ptr<portcullis::Junction> portcullis::Junction::parse(const string& line) {

    vector<string> parts; // #2: Search for tokens
    boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

    if (parts.size() != 36) {
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
            "Could not parse line due to incorrect number of columns. Expected 36 columns: ") + line));
    }

    // Create intron
    IntronPtr i(new Intron(
        RefSeq(lexical_cast<int32_t>(parts[1]), parts[2], lexical_cast<int32_t>(parts[3])),
        lexical_cast<int32_t>(parts[4]),
        lexical_cast<int32_t>(parts[5]),
        strandFromChar(parts[6][0])
    ));

    // Create basic junction
    shared_ptr<Junction> j(new Junction(
        i,
        lexical_cast<int32_t>(parts[7]),
        lexical_cast<int32_t>(parts[8])
    ));

    // Splice site strings
    j->setDa1(parts[9]);
    j->setDa2(parts[10]);

    // Set predictions to junction
    j->setPredictedStrand(strandFromChar(parts[11][0]));

    // Add metrics to junction
    j->setDonorAndAcceptorMotif(cssFromChar(parts[12][0]));
    j->setNbJunctionAlignments(lexical_cast<uint32_t>(parts[13]));
    j->setNbDistinctAlignments(lexical_cast<uint32_t>(parts[14]));
    j->setNbReliableAlignments(lexical_cast<uint32_t>(parts[15]));
    // Intron size not required
    j->setLeftAncSize(lexical_cast<uint32_t>(parts[17]));
    j->setRightAncSize(lexical_cast<uint32_t>(parts[18]));
    j->setMaxMinAnchor(lexical_cast<int32_t>(parts[19]));
    j->setDiffAnchor(lexical_cast<int32_t>(parts[20]));
    j->setNbDistinctAnchors(lexical_cast<uint32_t>(parts[21]));
    j->setEntropy(lexical_cast<double>(parts[22]));
    j->setMaxMMES(lexical_cast<uint32_t>(parts[23]));
    j->setHammingDistance5p(lexical_cast<uint32_t>(parts[24]));
    j->setHammingDistance3p(lexical_cast<uint32_t>(parts[25]));
    j->setCoverage(lexical_cast<double>(parts[26]));
    j->setUniqueJunction(lexical_cast<bool>(parts[27]));
    j->setPrimaryJunction(lexical_cast<bool>(parts[28]));
    j->setMultipleMappingScore(lexical_cast<double>(parts[29]));
    j->setNbMismatches(lexical_cast<uint16_t>(parts[30]));
    j->setNbMultipleSplicedReads(lexical_cast<uint32_t>(parts[31]));
    j->setNbUpstreamJunctions(lexical_cast<uint16_t>(parts[32]));
    j->setNbDownstreamJunctions(lexical_cast<uint16_t>(parts[33]));
    j->setNbUpstreamFlankingAlignments(lexical_cast<uint32_t>(parts[34]));
    j->setNbDownstreamFlankingAlignments(lexical_cast<uint32_t>(parts[35]));

    return j;
}
