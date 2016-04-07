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

#include <cstdint>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <unordered_map>
using std::boolalpha;
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

#include <portcullis/bam/bam_alignment.hpp>
#include <portcullis/bam/genome_mapper.hpp>
using portcullis::bam::CigarOp;
using portcullis::bam::GenomeMapper;
using portcullis::bam::Strand;

#include <portcullis/intron.hpp>
#include <portcullis/seq_utils.hpp>
using portcullis::Intron;

#include <portcullis/junction.hpp>


void portcullis::AlignmentInfo::calcMatchStats(const Intron& i, const uint32_t leftStart, const uint32_t rightEnd, const string& ancLeft, const string& ancRight) {
        
    uint32_t leftEnd = i.start - 1;
    uint32_t rightStart = i.end + 1;

    uint32_t qLeftStart = leftStart;
    uint32_t qLeftEnd = leftEnd;
    uint32_t qRightStart = rightStart;
    uint32_t qRightEnd = rightEnd;

    string qAnchorLeft = ba->getPaddedQuerySeq(leftStart, leftEnd, qLeftStart, qLeftEnd, false);
    string qAnchorRight = ba->getPaddedQuerySeq(rightStart, rightEnd, qRightStart, qRightEnd, false);

    string gAnchorLeft = ba->getPaddedGenomeSeq(ancLeft, leftStart, leftEnd, qLeftStart, qLeftEnd, false);
    string gAnchorRight = ba->getPaddedGenomeSeq(ancRight, rightStart, rightEnd, qRightStart, qRightEnd, false);

    if (qAnchorLeft.size() != gAnchorLeft.size() || qAnchorLeft.empty()) {
       BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
            "Left anchor region for query and genome are not the same size.")
               + "\nIntron: " + i.toString() 
               + "\nJunction anchor limits: " + lexical_cast<string>(leftStart) + "," + lexical_cast<string>(rightEnd)
               + "\nGenomic sequence: " + ancLeft
               + "\nAlignment coords (before soft clipping): " + ba->toString(false)
               + "\nAlignment coords (after soft clipping): " + ba->toString(true)
               + "\nRead name: " + ba->deriveName()
               + "\nRead seq (before soft clipping): " + ba->getQuerySeq() + " (" + lexical_cast<string>(ba->getQuerySeq().size()) + ")"
               + "\nRead seq (after soft clipping): " + ba->getQuerySeqAfterClipping() + " (" + lexical_cast<string>(ba->getQuerySeqAfterClipping().size()) + ")"
               + "\nCigar: " + ba->getCigarAsString()
               + "\nLeft Anchor query seq:  \n" + qAnchorLeft + " (" + lexical_cast<string>(qAnchorLeft.size()) + ")"
               + "\nLeft Anchor genome seq: \n" + gAnchorLeft + " (" + lexical_cast<string>(gAnchorLeft.size()) + ")"));
    }

    if (qAnchorRight.size() != gAnchorRight.size() || qAnchorRight.empty()) {
       BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
            "Right Anchor region for query and genome are not the same size.")
               + "\nIntron: " + i.toString() 
               + "\nJunction anchor limits: " + lexical_cast<string>(leftStart) + "," + lexical_cast<string>(rightEnd)
               + "\nGenomic sequence: " + ancRight
               + "\nAlignment coords (before soft clipping): " + ba->toString(false)
               + "\nAlignment coords (after soft clipping): " + ba->toString(true)
               + "\nRead name: " + ba->deriveName()
               + "\nRead seq (before soft clipping): " + ba->getQuerySeq() + " (" + lexical_cast<string>(ba->getQuerySeq().size()) + ")"
               + "\nRead seq (after soft clipping): " + ba->getQuerySeqAfterClipping() + " (" + lexical_cast<string>(ba->getQuerySeqAfterClipping().size()) + ")"
               + "\nCigar: " + ba->getCigarAsString()
               + "\nRight Anchor query seq:  \n" + qAnchorRight + " (" + lexical_cast<string>(qAnchorRight.size()) + ")"
               + "\nRight Anchor genome seq: \n" + gAnchorRight + " (" + lexical_cast<string>(gAnchorRight.size()) + ")"));
    }

    totalUpstreamMismatches = SeqUtils::hammingDistance(qAnchorLeft, gAnchorLeft);
    totalDownstreamMismatches = SeqUtils::hammingDistance(qAnchorRight, gAnchorRight);
    totalUpstreamMatches = qAnchorLeft.size() - totalUpstreamMismatches;
    totalDownstreamMatches = qAnchorRight.size() - totalDownstreamMismatches;
    nbMismatches = totalUpstreamMismatches + totalDownstreamMismatches;                       

    upstreamMatches = getNbMatchesFromEnd(qAnchorLeft, gAnchorLeft);
    downstreamMatches = getNbMatchesFromStart(qAnchorRight, gAnchorRight);
    minMatch = min(upstreamMatches, downstreamMatches);
    maxMatch = max(upstreamMatches, downstreamMatches);

    mmes = min(totalUpstreamMatches, totalDownstreamMatches);
}

uint32_t portcullis::AlignmentInfo::getNbMatchesFromStart(const string& query, const string& anchor) {

    for(size_t i = 0; i < query.size(); i++) {
        if (query[i] != anchor[i]) {
            return i;
        }
    }

    return query.size();
}

uint32_t portcullis::AlignmentInfo::getNbMatchesFromEnd(const string& query, const string& anchor) {
    for(size_t i = query.size() - 1; i >= 0; i--) {
        if (query[i] != anchor[i]) {
            return query.size() - i - 1;
        }
    }

    return query.size();
}

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
        return Strand::POSITIVE;
    }
    else if (seq == CANONICAL_SEQ_RC) {
        return Strand::NEGATIVE;
    }
    // Check for these after canonicals for performance reasons
    else if (seq == SEMI_CANONICAL_SEQ1 || seq == SEMI_CANONICAL_SEQ2) {
        return Strand::POSITIVE;
    }
    else if (seq == SEMI_CANONICAL_SEQ1_RC || seq == SEMI_CANONICAL_SEQ2_RC) {
        return Strand::NEGATIVE;
    }
    else {
        return Strand::UNKNOWN;
    }
}
    

portcullis::Junction::Junction(shared_ptr<Intron> _location, int32_t _leftAncStart, int32_t _rightAncEnd) :
    intron(_location) {
    leftAncStart = _leftAncStart;
    rightAncEnd = _rightAncEnd;
    meanQueryLength = 0;
    suspicious = false;
    pfp = false;
    canonicalSpliceSites = NO;
    maxMinAnchor = intron->minAnchorLength(_leftAncStart, _rightAncEnd);
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
    meanMismatches = 0;
    nbMultipleSplicedReads = 0;
    distanceToNextUpstreamJunction = 0;
    distanceToNextDownstreamJunction = 0;
    distanceToNearestJunction = 0;
    
    junctionOverhangs.clear();
    for(size_t i = 0; i < JO_NAMES.size(); i++) {
        junctionOverhangs.push_back(0);
    }
    
    readStrand = Strand::UNKNOWN;
    ssStrand = Strand::UNKNOWN;
    consensusStrand = Strand::UNKNOWN;
    
    nbUpstreamJunctions = 0;
    nbDownstreamJunctions = 0;
    
    alignments.clear();
    alignmentCodes.clear();
    trimmedCoverage.clear();
    trimmedLogDevCov.clear();
}
   
/**
 * Copy constructor, with option to copy bam alignments associated with the junction
 * @param j The other junction to deep copy into this
 * @param withAlignments Whether to copy over the alignments or not
 */
portcullis::Junction::Junction(const Junction& j, bool withAlignments) {
    intron = make_shared<Intron>(*(j.intron));        
    leftAncStart = j.leftAncStart;
    rightAncEnd = j.rightAncEnd;
    meanQueryLength = j.meanQueryLength;
    suspicious = j.suspicious;
    pfp = j.pfp;
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
    meanMismatches = j.meanMismatches;
    nbMultipleSplicedReads = j.nbMultipleSplicedReads;

    readStrand = j.readStrand;
    ssStrand = j.ssStrand;
    consensusStrand = j.consensusStrand;
    
    nbUpstreamJunctions = j.nbUpstreamJunctions;
    nbDownstreamJunctions = j.nbDownstreamJunctions;
    distanceToNextUpstreamJunction = j.distanceToNextUpstreamJunction;
    distanceToNextDownstreamJunction = j.distanceToNextDownstreamJunction;
    distanceToNearestJunction = j.distanceToNearestJunction;

    if (withAlignments) {
        for(size_t i = 0; i < j.alignments.size(); i++) {
            this->alignments.push_back(make_shared<AlignmentInfo>(j.alignments[i]->ba));
            this->alignmentCodes.push_back(j.alignments[i]->nameCode);
        }
    }
    
    trimmedCoverage.clear();
    for (auto& x : j.trimmedCoverage) {
        trimmedCoverage.push_back(x);
    }
    
    trimmedLogDevCov.clear();
    for (auto& x : j.trimmedLogDevCov) {
        trimmedLogDevCov.push_back(x);
    }
    
    junctionOverhangs.clear();
    for(size_t i = 0; i < JO_NAMES.size(); i++) {
        junctionOverhangs.push_back(j.getJunctionOverhangs(i));
    }
}
    
// **** Destructor ****
portcullis::Junction::~Junction() {
    alignments.clear();
    junctionOverhangs.clear();
}
    
void portcullis::Junction::clearAlignments() {
    alignments.clear();  
}
       
    
void portcullis::Junction::addJunctionAlignment(const BamAlignment& al) {

    // Make sure we take a proper copy of this alignment for safe storage
    AlignmentInfoPtr aip = make_shared<AlignmentInfo>(make_shared<BamAlignment>(al));
    this->alignments.push_back(aip);
    this->alignmentCodes.push_back(aip->nameCode);
    
    this->nbJunctionAlignments = this->alignments.size();

    if (al.getNbJunctionsInRead() > 1) {
        this->nbMultipleSplicedReads++;
    }
}
    
    
portcullis::CanonicalSS portcullis::Junction::setDonorAndAcceptorMotif(string seq1, string seq2) {
    this->canonicalSpliceSites = hasCanonicalSpliceSites(seq1, seq2);
    this->ssStrand = predictedStrandFromSpliceSites(seq1, seq2);
    // Also set consensusStrand here as readStrand should already be calculated
    this->consensusStrand = 
            this->readStrand == this->ssStrand ? this->readStrand : 
                this->readStrand == Strand::UNKNOWN ? this->ssStrand :
                    this->ssStrand == Strand::UNKNOWN ? this->readStrand :
                        Strand::UNKNOWN;
    
    this->da1 = this->consensusStrand == Strand::NEGATIVE ? SeqUtils::reverseComplement(seq1) : seq1;
    this->da2 = this->consensusStrand == Strand::NEGATIVE ? SeqUtils::reverseComplement(seq2) : seq2;
    
    return this->canonicalSpliceSites;
}
    
void portcullis::Junction::setFlankingAlignmentCounts(uint32_t nbUpstream, uint32_t nbDownstream) {
    this->nbUpstreamFlankingAlignments = nbUpstream;
    this->nbDownstreamFlankingAlignments = nbDownstream;
}
    

    
/**
 * Extends the anchor regions of the junction.  Also updates any relevant 
 * metrics.
 * @param otherStart The alternative start position of the left anchor
 * @param otherEnd The alternative end position of the right anchor
 */
void portcullis::Junction::extendAnchors(int32_t otherStart, int32_t otherEnd) {        

    leftAncStart = min(leftAncStart, otherStart);
    rightAncEnd = max(rightAncEnd, otherEnd);

    int32_t otherMinAnchor = intron->minAnchorLength(otherStart, otherEnd);

    maxMinAnchor = max(maxMinAnchor, otherMinAnchor);
}
        

void portcullis::Junction::determineStrandFromReads() {
    
    uint32_t nb_pos = 0;
    uint32_t nb_neg = 0;
    uint32_t nb_unk = 0;
    for(const auto& a : alignments) {
        switch (a->ba->getStrand()) {
            case Strand::POSITIVE:
                nb_pos++;
                break;
            case Strand::NEGATIVE:
                nb_neg++;
                break;
            case Strand::UNKNOWN:
                nb_unk++;
                break;
        }        
    }
    
    uint32_t total = nb_pos + nb_neg + nb_unk;
    const double threshold = 0.95;
    
    if ((double)nb_pos / (double)total >= threshold) {
        readStrand = Strand::POSITIVE;
    }
    else if ((double)nb_neg / (double)total >= threshold) {
        readStrand = Strand::NEGATIVE;
    }
    else {
        readStrand = Strand::UNKNOWN;
    }
}
    
void portcullis::Junction::processJunctionWindow(const GenomeMapper& genomeMapper) {

    if (intron == nullptr) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find genomic sequence for this junction as no intron is defined")));
    
    // Process the predicted donor / acceptor regions and update junction
    int donorLen = -1;
    int acceptorLen = -1;
    string donor = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->start, intron->start + 1, &donorLen);
    string acceptor = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->end - 1, intron->end, &acceptorLen);
    
    if (donorLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find donor site (left side splice site) region for junction: ") + this->intron->toString()));
    
    if (donorLen != 2)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence for left side splice site of junction ") + this->intron->toString() + " is not the expected length" +
                "\nRetrieved sequence Length: " + lexical_cast<string>(donorLen) + 
                "\nExpected sequence length: " + lexical_cast<string>(2) + "\n"));
        
    if (acceptorLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find acceptor site (right side splice site) region for junction: ") + this->intron->toString()));
    
    if (acceptorLen != 2)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence for right side splice site of junction ") + this->intron->toString() + " is not the expected length" +
                "\nRetrieved sequence Length: " + lexical_cast<string>(acceptorLen) + 
                "\nExpected sequence length: " + lexical_cast<string>(2) + "\n"));
        
    boost::to_upper(donor);    // Removes any lowercase bases representing repeats
    boost::to_upper(acceptor);    // Removes any lowercase bases representing repeats
    this->setDonorAndAcceptorMotif(donor, acceptor);
    
    // Just access the whole junction region
    int leftAncLen = -1;
    int leftIntLen = -1;
    int rightAncLen = -1;
    int rightIntLen = -1;
    string leftAnc = genomeMapper.fetchBases(intron->ref.name.c_str(), leftAncStart, intron->start - 1, &leftAncLen);
    string rightAnc = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->end + 1, rightAncEnd, &rightAncLen);
    string leftInt = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->start, intron->start + 9, &leftIntLen);
    string rightInt = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->end - 9, intron->end, &rightIntLen);

    if (leftAncLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find left anchor region for junction: ") + this->intron->toString()));

    if (rightAncLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find right anchor region for junction: ") + this->intron->toString()));
    
    if (leftIntLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find left intron region for junction: ") + this->intron->toString()));

    if (rightIntLen == -1) 
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Can't find right intron region for junction: ") + this->intron->toString()));

    int expLeftLen = intron->start - leftAncStart;
    if (leftAncLen != expLeftLen && expLeftLen > 0)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence for left anchor of junction ") + this->intron->toString() + " is not the expected length" +
                "\nRetrieved sequence Length: " + lexical_cast<string>(leftAncLen) + 
                "\nExpected sequence length: " + lexical_cast<string>(expLeftLen) + 
                "\nRetrieved sequence: " + leftAnc));

    int expRightLen = rightAncEnd - intron->end;
    if (rightAncLen != expRightLen && expRightLen > 0)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence for right anchor of junction ") + this->intron->toString() + " is not the expected length" +                
                "\nRetrieved sequence Length: " + lexical_cast<string>(rightAncLen) + 
                "\nExpected sequence length: " + lexical_cast<string>(expRightLen) + 
                "\nRetrieved sequence: " + rightAnc));
        
    if (leftIntLen != 10)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence for left intron region of junction ") + this->intron->toString() + " is not the expected length" +
                "\nRetrieved sequence Length: " + lexical_cast<string>(leftIntLen) + 
                "\nExpected sequence length: " + lexical_cast<string>(10) + "\n"));

    if (rightIntLen != 10)
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Retrieved sequence for right intron region of junction ") + this->intron->toString() + " is not the expected length" +
                "\nRetrieved sequence Length: " + lexical_cast<string>(rightIntLen) + 
                "\nExpected sequence length: " + lexical_cast<string>(10) + "\n"));

    // Removes any lowercase bases representing repeats
    boost::to_upper(leftAnc);
    boost::to_upper(rightAnc);
    boost::to_upper(leftInt);
    boost::to_upper(rightInt);
    
    string leftAnchor10 = leftAncLen < 10 ? leftAnc : leftAnc.substr(leftAncLen - 10, 10);
    string rightAnchor10 = rightAncLen < 10 ? rightAnc : rightAnc.substr(0, 10);
    
    this->calcHammingScores(leftAnchor10, leftInt, rightInt, rightAnchor10);
    
    // Update match statistics for each alignment
    for(const auto& a : alignments) {
        a->calcMatchStats(*getIntron(), this->getLeftAncStart(), this->getRightAncEnd(), leftAnc, rightAnc);
    }
    
    // MaxMMES can now use info in alignments
    this->calcMismatchStats();    
}
    
void portcullis::Junction::processJunctionVicinity(BamReader& reader, int32_t refLength, int32_t meanQueryLength, int32_t maxQueryLength) {

    int32_t refId = intron->ref.index;
    
    uint32_t nbLeftFlankingAlignments = 0, nbRightFlankingAlignments = 0;

    int32_t regionStart = leftAncStart - maxQueryLength - 1; regionStart < 0 ? 0 : regionStart;
    int32_t regionEnd = rightAncEnd + maxQueryLength + 1; regionEnd >= refLength ? refLength - 1 : regionEnd;

    // Focus only on the (expanded... to be safe...) region of interest
    reader.setRegion(refId, regionStart, regionEnd);

    while(reader.next()) {

        const BamAlignment& ba = reader.current();
        
        int32_t pos = ba.getStart();

        //TODO: Should we consider strand specific reads differently here?

        // Look for left flanking alignments
        if (    intron->start > pos && 
                leftAncStart <= ba.getEnd()) {
            nbLeftFlankingAlignments++;
        }

        // Look for right flanking alignments
        if (    rightAncEnd >= pos && 
                intron->end < pos) {
            nbRightFlankingAlignments++;
        }
    }

    this->setFlankingAlignmentCounts(nbLeftFlankingAlignments, nbRightFlankingAlignments);
}
    
    
void portcullis::Junction::calcMetrics() {
    
    determineStrandFromReads();
    calcAnchorStats();      // Metrics 5 and 7
    calcEntropy();          // Metric 6
    calcAlignmentStats();   // Metrics 8, 9 and 19
}


    
/**
 * Metric 5 and 7: Diff Anchor and # Distinct Anchors
 * @return 
 */
void portcullis::Junction::calcAnchorStats() {

    int32_t minLeftSize = INT32_MAX, minRightSize = INT32_MAX; 
    int32_t maxLeftSize = 0, maxRightSize = 0;
    int32_t lastLStart = -1, lastREnd = -1;

    nbDistinctAnchors = 0;    
    
    for(const auto& a : alignments) {
        
        BamAlignmentPtr ba = a->ba;

        const int32_t lStart = ba->getStart();
        const int32_t rEnd = ba->getEnd();
        
        int32_t leftSize = ba->calcNbAlignedBases(lStart, intron->start - 1, false);
        int32_t rightSize = ba->calcNbAlignedBases(intron->end + 1, rEnd, false);   // Will auto crop for end of alignment

        maxLeftSize = max(maxLeftSize, leftSize);
        minLeftSize = min(minLeftSize, leftSize);
        maxRightSize = max(maxRightSize, rightSize);
        minRightSize = min(minRightSize, rightSize);
        
        if (lStart != lastLStart || rEnd != lastREnd) {
            nbDistinctAnchors++;
            lastLStart = lStart;
            lastREnd = rEnd;
        }        
    }

    int32_t diffLeftSize = maxLeftSize - minLeftSize; diffLeftSize = diffLeftSize < 0 ? 0 : diffLeftSize;
    int32_t diffRightSize = maxRightSize - minRightSize; diffLeftSize = diffRightSize < 0 ? 0 : diffRightSize;   
    
    leftAncSize = maxLeftSize;
    rightAncSize = maxRightSize;            
    diffAnchor = min(diffLeftSize, diffRightSize);    
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
double portcullis::Junction::calcEntropy() {

    vector<int32_t> junctionPositions;

    for(const auto& a : alignments) {
        junctionPositions.push_back(a->ba->getStart());
    }
    
    // Should already be sorted but let's be sure.  This is critical to the rest
    // of the algorithm.  It's possible after soft clips are removed that the reads
    // are not strictly in the correct order.
    std::sort(junctionPositions.begin(), junctionPositions.end());
    
    return calcEntropy(junctionPositions);
}
    
double portcullis::Junction::calcEntropy(const vector<int32_t> junctionPositions) {
    
    size_t nbJunctionAlignments = junctionPositions.size();

    if (nbJunctionAlignments <= 1)
        return 0;

    double sum = 0.0;

    int32_t lastOffset = junctionPositions[0];
    uint32_t readsAtOffset = 0;

    for(size_t i = 0; i < nbJunctionAlignments; i++) {

        int32_t pos = junctionPositions[i];

        readsAtOffset++;

        if (pos != lastOffset || i == nbJunctionAlignments - 1) {
            double pI = (double)readsAtOffset / (double)nbJunctionAlignments;
            sum += pI * log2(pI);
            lastOffset = pos;
            readsAtOffset = 0;
        }
    }

    entropy = fabs(sum);
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
    nbUpstreamJunctions = 0;
    nbDownstreamJunctions = 0;

    //cout << junctionAlignments.size() << endl;

    for(const auto& a : alignments) {

        BamAlignmentPtr ba = a->ba;
        
        const int32_t start = ba->getStart();
        const int32_t end = ba->getEnd();

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
                else if (pos > intron->end + 1) {
                    downjuncs++;
                }  
            }
        }

        nbUpstreamJunctions = max(nbUpstreamJunctions, upjuncs);
        nbDownstreamJunctions = max(nbDownstreamJunctions, downjuncs);
    }        
}
    
/**
 * Metric 13 and 14: Calculates the 5' and 3' hamming distances from a genomic
 * region represented by this junction
 * 
 * Note, the logic in here is quite complex (and initially looks weird) but it is 
 * correct... only change it if you know what you are doing!
 */
void portcullis::Junction::calcHammingScores(   const string& leftAnchor, const string& leftIntron, 
                                                const string& rightIntron, const string& rightAnchor) {

    const int32_t leftDelta = leftAnchor.size() - rightIntron.size();
    const int32_t leftOffset = leftDelta <= 0 ? 0 : leftDelta;
    const int32_t leftLen = min(leftAnchor.size(), rightIntron.size());
    const int32_t rightLen = min(leftIntron.size(), rightAnchor.size());
    
    // TODO, might want to modify this logic later, but worst case is that the 5' and
    // 3' results are swapped
    Strand s = consensusStrand != UNKNOWN ? consensusStrand : Strand::UNKNOWN;

    const string la = leftAnchor.size() > leftLen ? leftAnchor.substr(leftOffset, leftLen) : leftAnchor;
    const string li = leftIntron.size() > rightLen ? leftIntron.substr(0, rightLen) : leftIntron;
    const string ri = rightIntron.size() > leftLen ? rightIntron.substr(leftOffset, leftLen) : rightIntron;
    const string ra = rightAnchor.size() > rightLen ? rightAnchor.substr(0, rightLen) : rightAnchor;

    string anchor5p;
    string intron5p;
    string intron3p;
    string anchor3p;
    
    if (s == NEGATIVE) {
        anchor5p = SeqUtils::reverseComplement(ra);
        intron5p = SeqUtils::reverseComplement(ri);
        intron3p = SeqUtils::reverseComplement(li);
        anchor3p = SeqUtils::reverseComplement(la);
    }
    else {
        // This is for the positive strand or if strand is unknown (i.e. novel splice site)
        // This shouldn't make much difference for downstream analysis.  Worst case is that
        // hamming 5' and hamming 3' figures are reversed.
        anchor5p = la;
        intron5p = li;
        intron3p = ri;
        anchor3p = ra;
    }
     
    hammingDistance5p = SeqUtils::hammingDistance(anchor5p, intron3p);
    hammingDistance3p = SeqUtils::hammingDistance(anchor3p, intron5p);  
}
    
/**
 * Calculates MaxMMES, mismatches, and junction overhangs.
 */
void portcullis::Junction::calcMismatchStats() {

    uint32_t nbMismatches = 0;
    uint32_t firstMismatch = 100000000;
    for(const auto& a : alignments) {
        
        // Update maxMMES for this alignment
        maxMMES = max(maxMMES, a->mmes);
        
        // Update total number of mismatches in this junction
        nbMismatches += a->nbMismatches;
        
        // Keep a record of the first mismatch detected
        if (a->minMatch > 0) {
            firstMismatch = min(firstMismatch, a->minMatch);
        }
        
        // Update junction overhang vector
        for(uint16_t i = 0; i < JO_NAMES.size() && i < a->minMatch; i++) {
            junctionOverhangs[i]++;
        }
    }
    
    // Set mean mismatches across junction
    meanMismatches = (double)nbMismatches / (double)alignments.size();
    
    // Assuming we have some mismatches determine if this junction has no overhangs
    // extending beyond first mismatch.  If so determine if that distance is small
    // enough to consider the junction as suspicious
    if (nbMismatches > 0) {
        bool found = false;
        for(const auto& a : alignments) {
            if (a->mmes > firstMismatch) {
                found = true;
                break;
            }
        }
        
        if (!found) {
            suspicious = true;
        }
    }
}

    
/**
 * Calculates metric 18.  Multiple mapping score
 */
void portcullis::Junction::calcMultipleMappingScore(SplicedAlignmentMap& map) {

    size_t N = alignmentCodes.size();

    uint32_t M = 0;
    for(const auto& a : alignmentCodes) {
        M += map[a];  // Number of multiple splitting patterns
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
         << "Anchor limits: (" << leftAncStart << ", " << rightAncEnd << ")" << delimiter
         << "Junction Predictions:" << delimiter
         << "Reads Strand: " << strandToString(readStrand) << delimiter
         << "Splice Site Strand: " << strandToString(ssStrand) << delimiter
         << "Consensus Strand: " << strandToString(consensusStrand) << delimiter
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
         << "M19: Mean mismatches: " << meanMismatches << delimiter
         << "M20: # Multiple Spliced Reads: " << nbMultipleSplicedReads << delimiter
         << "M21: # Upstream Junctions: " << nbUpstreamJunctions << delimiter
         << "M22: # Downstream Junctions: " << nbDownstreamJunctions << delimiter
         << "M23: # Upstream Non-Spliced Alignments: " << nbUpstreamFlankingAlignments << delimiter
         << "M24: # Downstream Non-Spliced Alignments: " << nbDownstreamFlankingAlignments << delimiter
         << "M25: Distance to next upstream junction: " << distanceToNextUpstreamJunction << delimiter
         << "M26: Distance to next downstream junction: " << distanceToNextDownstreamJunction << delimiter
         << "M27: Distance to nearest junction: " << distanceToNearestJunction;         
}
    
/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::condensedOutputDescription(std::ostream &strm, string delimiter) {

    strm << "Reads Strand: " << strandToString(readStrand) << delimiter
         << "Splice Site Strand: " << strandToString(ssStrand) << delimiter
         << "Consensus Strand: " << strandToString(consensusStrand) << delimiter
         << "M1-Canonical?=" << cssToString(canonicalSpliceSites) << delimiter
         << "M2-NbAlignments=" << getNbJunctionAlignments() << delimiter
         << "M3-NbDistinctAlignments=" << nbDistinctAlignments << delimiter
         << "M4-NbReliableAlignments=" << nbReliableAlignments << delimiter
         << "M5-IntronSize=" << getIntronSize() << delimiter
         << "M6-LeftAnchorSize=" << leftAncSize << delimiter
         << "M7-RightAnchorSize=" << rightAncSize << delimiter
         << "M8-MaxMinAnchor=" << maxMinAnchor << delimiter
         << "M9-DiffAnchor=" << diffAnchor << delimiter
         << "M10-NbDistinctAnchors=" << nbDistinctAnchors << delimiter
         << "M11-Entropy=" << entropy << delimiter
         << "M12-MaxMMES=" << maxMMES << delimiter
         << "M13-HammingDistance5=" << hammingDistance5p << delimiter
         << "M14-HammingDistance3=" << hammingDistance3p << delimiter
         << "M15-Coverage=" << coverage << delimiter
         << "M16-UniqueJunction=" << boolalpha << uniqueJunction << delimiter
         << "M17-PrimaryJunction=" << boolalpha << primaryJunction << delimiter
         << "M18-MultipleMappingScore=" << multipleMappingScore << delimiter
         << "M19-MeanMismatches=" << meanMismatches << delimiter
         << "M20-NbMultipleSplicedReads=" << nbMultipleSplicedReads << delimiter
         << "M21-NbUpstreamJunctions=" << nbUpstreamJunctions << delimiter
         << "M22-NbDownstreamJunctions=" << nbDownstreamJunctions << delimiter
         << "M23-NbUpstreamNonSplicedAlignments=" << nbUpstreamFlankingAlignments << delimiter
         << "M24-NbDownstreamNonSplicedAlignments=" << nbDownstreamFlankingAlignments << delimiter
         << "M25-DistanceToNextUpstreamJunction=" << distanceToNextUpstreamJunction << delimiter
         << "M26-DistanceToNextDownstreamJunction=" << distanceToNextDownstreamJunction << delimiter
         << "M27-DistanceToNearestJunction=" << distanceToNearestJunction << delimiter
         << "PFP=" << boolalpha << pfp;         
}

    
/**
 * Complete human readable description of this intron (for augustus hints)
 * @param strm
 */
void portcullis::Junction::outputIntronGFF(std::ostream &strm, uint32_t id, const string& source) {

    // Use intron strand if known, otherwise use the predicted strand,
    // if predicted strand is also unknown then use "." to indicated unstranded
    const char strand = consensusStrand == UNKNOWN ? '.' : strandToChar(consensusStrand);

    string juncId = string("junc_") + lexical_cast<string>(id);

    // Modify coordinates to 1-based end inclusive
    
    // Output junction parent
    strm << intron->ref.name << "\t"
         << source << "\t"    // source
         << "intron" << "\t"        // type (may change later)
         << intron->start + 1 << "\t"   // start
         << intron->end  + 1<< "\t"     // end
         << nbJunctionAlignments << "\t"           // No score for the moment
         << strand << "\t"          // strand
         << "." << "\t"             // Just put "." for the phase
         // Removing this as it causes issues with PASA downstream
         /**<< "Note=cov:" << nbJunctionAlignments 
                        << "|rel:" << this->nbReliableAlignments 
                        << "|ent:" << std::setprecision(4) << this->entropy << std::setprecision(9) 
                        << "|maxmmes:" << this->maxMMES
                        << "|ham:" << min(this->hammingDistance3p, this->hammingDistance5p) << ";"  // Number of times it was seen**/
         << "mult=" << nbJunctionAlignments << ";"  // Coverage for augustus
         << "grp=" << juncId << ";"  // ID for augustus
         << "src=E";                // Source for augustus
    strm << endl;

}

/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::outputJunctionGFF(std::ostream &strm, uint32_t id, const string& source) {

    // Use intron strand if known, otherwise use the predicted strand,
    // if predicted strand is also unknown then use "." to indicated unstranded
    const char strand = consensusStrand == Strand::UNKNOWN ? '.' : strandToChar(consensusStrand);

    string juncId = string("junc_") + lexical_cast<string>(id);

    // Modify coordinates to 1-based end inclusive
    
    // Output junction parent
    strm << intron->ref.name << "\t"
         << source << "\t"    // source
         << "match" << "\t"      // type (may change later)
         << leftAncStart + 1 << "\t"  // start
         << rightAncEnd + 1 << "\t"   // end
         << "0.0" << "\t"           // No score for the moment
         << strand << "\t"          // strand
         << "." << "\t"             // Just put "." for the phase
         << "ID=" << juncId << ";"  // ID of the intron
         << "Note=cov:" << nbJunctionAlignments 
                        << "|rel:" << this->nbReliableAlignments 
                        << "|ent:" << std::setprecision(4) << this->entropy << std::setprecision(9)
                        << "|maxmmes:" << this->maxMMES
                        << "|ham:" << min(this->hammingDistance3p, this->hammingDistance5p) << ";"  // Number of times it was seen
         << "mult=" << nbJunctionAlignments << ";"  // Coverage for augustus
         << "grp=" << juncId << ";"  // ID for augustus
         << "src=E;";                // Source for augustus
    condensedOutputDescription(strm, ";");
    strm << endl;

    // Make modifications to coordinates so they are suitable for GFF.  1-based with end positions inclusive.
    
    // Output left exonic region
    strm << intron->ref.name << "\t"
         << source << "\t"
         << "match_part" << "\t"
         << leftAncStart + 1 << "\t"
         << (intron->start) << "\t"
         << "0.0" << "\t"
         << strand << "\t"
         << "." << "\t"
         << "ID=" << juncId << "_left" << ";"
         << "Parent=" << juncId << endl;

    // Output right exonic region
    strm << intron->ref.name << "\t"
         << source << "\t"
         << "match_part" << "\t"
         << (intron->end + 2) << "\t"
         << rightAncEnd + 1 << "\t"
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
void portcullis::Junction::outputBED(std::ostream &strm, const string& prefix, uint32_t id) {

    // Use intron strand if known, otherwise use the predicted strand,
    // if predicted strand is also unknown then use "." to indicated unstranded
    const char strand = consensusStrand == Strand::UNKNOWN ?
                                '.' :
                                strandToChar(consensusStrand);

    string juncId = prefix + "_" + lexical_cast<string>(id);

    int32_t sz1 = intron->start - leftAncStart;
    int32_t sz2 = rightAncEnd - intron->end;
    string blockSizes = lexical_cast<string>(sz1) + "," + lexical_cast<string>(sz2);
    string blockStarts = lexical_cast<string>(0) + "," + lexical_cast<string>(intron->end - leftAncStart + 1);

    // Output junction parent
    strm << intron->ref.name << "\t"         // chrom
         << leftAncStart << "\t"  // chromstart
         << rightAncEnd + 1 << "\t"   // chromend (adding 1 as end position is exclusive)
         << juncId << "\t"          // name
         << this->getNbJunctionAlignments() << "\t"           // Use the depth as the score for the moment
         << strand << "\t"          // strand
         << intron->start << "\t"   // thickstart
         << intron->end + 1 << "\t"     // thickend  (adding 1 as end position is exclusive)
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
            boost::algorithm::join(STRAND_NAMES, "\t") + "\t" +
            boost::algorithm::join(METRIC_NAMES, "\t") + "\t" + 
            "MQL\tSuspect\tPFP\t" +
            boost::algorithm::join(JO_NAMES, "\t");
}

shared_ptr<portcullis::Junction> portcullis::Junction::parse(const string& line) {

    vector<string> parts; // #2: Search for tokens
    boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

    uint32_t expected_cols = 13 + STRAND_NAMES.size() + METRIC_NAMES.size() + JO_NAMES.size();
    
    if (parts.size() != expected_cols) {
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
            "Could not parse line due to incorrect number of columns.  This is probably a version mismatch.  Check file and portcullis versions.  Expected ") 
                + std::to_string(expected_cols) + " columns.  Found " 
                + std::to_string(parts.size()) + ".  Line: " + line));
    }

    // Create intron
    IntronPtr i = make_shared<Intron>(
        RefSeq(
            lexical_cast<int32_t>(parts[1]), 
            parts[2], 
            lexical_cast<int32_t>(parts[3])
        ),
        lexical_cast<int32_t>(parts[4]),
        lexical_cast<int32_t>(parts[5])
    );

    // Create basic junction
    shared_ptr<Junction> j = make_shared<Junction>(
        i,
        lexical_cast<int32_t>(parts[6]),
        lexical_cast<int32_t>(parts[7])
    );

    // Splice site strings
    j->setDa1(parts[8]);
    j->setDa2(parts[9]);

    // Set predictions to junction
    j->readStrand = strandFromChar(parts[10][0]);
    j->ssStrand = strandFromChar(parts[11][0]);
    j->consensusStrand = strandFromChar(parts[12][0]);

    // Add metrics to junction
    j->setDonorAndAcceptorMotif(cssFromChar(parts[13][0]));
    j->setNbJunctionAlignments(lexical_cast<uint32_t>(parts[14]));
    j->setNbDistinctAlignments(lexical_cast<uint32_t>(parts[15]));
    j->setNbReliableAlignments(lexical_cast<uint32_t>(parts[16]));
    // Intron size not required
    j->setLeftAncSize(lexical_cast<uint32_t>(parts[18]));
    j->setRightAncSize(lexical_cast<uint32_t>(parts[19]));
    j->setMaxMinAnchor(lexical_cast<int32_t>(parts[20]));
    j->setDiffAnchor(lexical_cast<int32_t>(parts[21]));
    j->setNbDistinctAnchors(lexical_cast<uint32_t>(parts[22]));
    j->setEntropy(lexical_cast<double>(parts[23]));
    j->setMaxMMES(lexical_cast<uint32_t>(parts[24]));
    j->setHammingDistance5p(lexical_cast<uint32_t>(parts[25]));
    j->setHammingDistance3p(lexical_cast<uint32_t>(parts[26]));
    j->setCoverage(lexical_cast<double>(parts[27]));
    j->setUniqueJunction(lexical_cast<bool>(parts[28]));
    j->setPrimaryJunction(lexical_cast<bool>(parts[29]));
    j->setMultipleMappingScore(lexical_cast<double>(parts[30]));
    j->setMeanMismatches(lexical_cast<double>(parts[31]));
    j->setNbMultipleSplicedReads(lexical_cast<uint32_t>(parts[32]));
    j->setNbUpstreamJunctions(lexical_cast<uint16_t>(parts[33]));
    j->setNbDownstreamJunctions(lexical_cast<uint16_t>(parts[34]));
    j->setNbUpstreamFlankingAlignments(lexical_cast<uint32_t>(parts[35]));
    j->setNbDownstreamFlankingAlignments(lexical_cast<uint32_t>(parts[36]));
    j->setDistanceToNextUpstreamJunction(lexical_cast<uint32_t>(parts[37]));
    j->setDistanceToNextDownstreamJunction(lexical_cast<uint32_t>(parts[38]));
    j->setDistanceToNearestJunction(lexical_cast<uint32_t>(parts[39]));
    
    // Read Junction overhangs
    j->setMeanQueryLength(lexical_cast<uint32_t>(parts[40]));
    j->setSuspicious(lexical_cast<bool>(parts[41]));
    j->setPotentialFalsePositive(lexical_cast<bool>(parts[42]));
    for(size_t i = 0; i < JO_NAMES.size(); i++) {
        j->setJunctionOverhangs(i, lexical_cast<uint32_t>(parts[43 + i]));
    }
    
    return j;
}
