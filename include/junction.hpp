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
using std::endl;
using std::min;
using std::max;
using std::string;
using std::size_t;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
using boost::lexical_cast;
using boost::shared_ptr;

#include <api/BamAlignment.h>
#include <utils/bamtools_pileup_engine.h>
using namespace BamTools;

#include "intron.hpp"
#include "genome_mapper.hpp"
#include "bam_utils.hpp"
#include "seq_utils.hpp"
using portculis::Intron;
using portculis::Strand;
using portculis::bamtools::BamUtils;

namespace portculis {    

const uint16_t MAP_QUALITY_THRESHOLD = 30;    
    

typedef boost::error_info<struct JunctionError,string> JunctionErrorInfo;
struct JunctionException: virtual boost::exception, virtual std::exception { };

const string METRIC_NAMES[18] = {
        "M1-nb_reads",
        "M2-canonical_ss",
        "M3-intron_size",
        "M4-max_min_anc",
        "M5-dif_anc",
        "M6-entropy",
        "M7-dist_anc",
        "M8-nb_dist_aln",
        "M9-nb_rel_aln",
        "M10-up_aln",
        "M11-down_aln",
        "M12-maxmmes",
        "M13-hamming5p",
        "M14-hamming3p",
        "M15-coverage",
        "M16-uniq_junc",
        "M17-primary_junc",
        "M18-mm_score"
    };

const string PREDICTION_NAMES[1] = {
        "P1-strand"
    };

const string CANONICAL_SEQ = "GTAG";
const string SEMI_CANONICAL_SEQ1 = "ATAC";
const string SEMI_CANONICAL_SEQ2 = "GCAG";

const string CANONICAL_SEQ_RC = SeqUtils::reverseComplement(CANONICAL_SEQ);
const string SEMI_CANONICAL_SEQ1_RC = SeqUtils::reverseComplement(SEMI_CANONICAL_SEQ1);
const string SEMI_CANONICAL_SEQ2_RC = SeqUtils::reverseComplement(SEMI_CANONICAL_SEQ2);

enum CanonicalSS {
    CANONICAL,
    SEMI_CANONICAL,
    NO
};

static CanonicalSS cssFromChar(char css) {
    switch(css) {
        case 'C':
            return CANONICAL;
        case 'S':
            return SEMI_CANONICAL;
        case 'N':
            return NO;
    }

    return NO;
}

static char cssToChar(CanonicalSS css) {
    
    switch(css) {
        case CANONICAL:
            return 'C';
        case SEMI_CANONICAL:
            return 'S';
        case NO:
            return 'N';
    }

    return 'N';
}

static string cssToString(CanonicalSS css) {
    
    switch(css) {
        case CANONICAL:
            return "Canonical";
        case SEMI_CANONICAL:
            return "Semi-canonical";
        case NO:
            return "No";
    }

    return string("No");
}

class Junction {
    
    
private:
    
    // **** Properties that describe where the junction is ****
    shared_ptr<Intron> intron;
    vector<BamAlignment> junctionAlignments;
    
    
    // **** Junction metrics ****
                                                // Metric 1 (nbReads) derived from size of junction alignment vector
    CanonicalSS canonicalSpliceSites;           // Metric 2
                                                // Metric 3 (intron size) calculated via location properties
    int32_t  maxMinAnchor;                      // Metric 4
    int32_t  diffAnchor;                        // Metric 5
    double   entropy;                           // Metric 6
    uint32_t nbDistinctAnchors;                 // Metric 7
    uint32_t nbDistinctAlignments;              // Metric 8
    uint32_t nbReliableAlignments;              // Metric 9
    uint32_t nbUpstreamFlankingAlignments;      // Metric 10
    uint32_t nbDownstreamFlankingAlignments;    // Metric 11
    uint32_t maxMMES;                           // Metric 12
    int16_t  hammingDistance5p;                 // Metric 13
    int16_t  hammingDistance3p;                 // Metric 14
    double   coverage;                          // Metric 15
    bool     uniqueJunction;                    // Metric 16
    bool     primaryJunction;                   // Metric 17
    double   multipleMappingScore;              // Metric 18
    
    
    // **** Predictions ****
    
    Strand predictedStrand;
    
    
    
    // **** Additional properties ****
    
    int32_t leftFlankStart;
    int32_t rightFlankEnd;
    string da1, da2;                    // These store the nucleotides found at the predicted donor / acceptor sites in the intron
    
    
    
    void init(int32_t _leftFlankStart, int32_t _rightFlankEnd) {
        
        leftFlankStart = _leftFlankStart;
        rightFlankEnd = _rightFlankEnd;
        canonicalSpliceSites = NO;
        maxMinAnchor = intron->minAnchorLength(_leftFlankStart, _rightFlankEnd);
        diffAnchor = 0;
        entropy = 0;
        nbDistinctAnchors = 0;
        nbDistinctAlignments = 0;
        nbReliableAlignments = 0;
        nbUpstreamFlankingAlignments = 0;
        nbDownstreamFlankingAlignments = 0;
        maxMMES = 0;
        hammingDistance5p = -1;
        hammingDistance3p = -1;
        coverage = 0.0;
        uniqueJunction = false;
        primaryJunction = false;
        multipleMappingScore = 0.0;
        
        predictedStrand = UNKNOWN;
    }
    
    
protected:
    
    /**
     * Tests whether the two strings could represent valid donor and acceptor sites
     * for this junction
     * @param seq1
     * @param seq2
     * @return 
     */
    CanonicalSS hasCanonicalSpliceSites(const string& seq1, const string& seq2) {
        
        if (intron == NULL || seq1.size() != 2 || seq2.size() != 2)
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
    
    Strand predictedStrandFromSpliceSites(const string& seq1, const string& seq2) {
        
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
    
    
public:
    
    // **** Constructors ****
    
    Junction(shared_ptr<Intron> _location, int32_t _leftFlankStart, int32_t _rightFlankEnd) :
        intron(_location) {
        init(_leftFlankStart, _rightFlankEnd);
    }
    
    // **** Destructor ****
    virtual ~Junction() {
        
    }
    
   
    shared_ptr<Intron> getLocation() const {
        return intron;
    }

    void setLocation(shared_ptr<Intron> location) {
        this->intron = location;
    }

    
    
    void addJunctionAlignment(const BamAlignment& al) {
        this->junctionAlignments.push_back(al);                        
    }
    
    void setDonorAndAcceptorMotif(CanonicalSS canonicalSS) {
        this->canonicalSpliceSites = canonicalSS;
    }
    
    CanonicalSS setDonorAndAcceptorMotif(string seq1, string seq2) {
        this->da1 = seq1;
        this->da2 = seq2;
        this->canonicalSpliceSites = hasCanonicalSpliceSites(seq1, seq2);
        this->predictedStrand = predictedStrandFromSpliceSites(seq1, seq2);
        return this->canonicalSpliceSites;
    }
    
    void setFlankingAlignmentCounts(uint32_t nbUpstream, uint32_t nbDownstream) {
        this->nbUpstreamFlankingAlignments = nbUpstream;
        this->nbDownstreamFlankingAlignments = nbDownstream;
    }
    
    void setPrimaryJunction(bool primaryJunction) {
        this->primaryJunction = primaryJunction;
    }

    void setUniqueJunction(bool uniqueJunction) {
        this->uniqueJunction = uniqueJunction;
    }

    
    bool sharesDonorOrAcceptor(shared_ptr<Junction> other) {
        return this->intron->sharesDonorOrAcceptor(*(other->intron));
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
        
        int32_t otherMinAnchor = intron->minAnchorLength(otherStart, otherEnd);
        
        maxMinAnchor = max(maxMinAnchor, otherMinAnchor);
    }
        
    
    int32_t processJunctionWindow(GenomeMapper* genomeMapper, RefVector& refs) {
        
        if (intron == NULL) 
            BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                    "Can't find genomic sequence for this junction as no intron is defined")));
                
        int32_t refid = intron->refId;
        char* refName = new char[refs[refid].RefName.size() + 1];
        strcpy(refName, refs[refid].RefName.c_str());

        // Just access the whole junction region
        int seqLen = -1;
        string region(genomeMapper->fetchBases(refName, leftFlankStart, rightFlankEnd, &seqLen));        
        if (seqLen == -1) 
            BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                    "Can't find genomic region for junction")));
        if (seqLen != rightFlankEnd - leftFlankStart + 1)
            BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                    "Retrieved sequence is not of the expected length")));
            

        // Process the predicted donor / acceptor regions and update junction
        string daSeq1 = region.substr(intron->start - leftFlankStart, 2);
        string daSeq2 = region.substr(intron->end - 2 - leftFlankStart, 2);
        bool validDA = setDonorAndAcceptorMotif(daSeq1, daSeq2);
           
        // Create strings for hamming analysis
        calcHammingScores(region);
        
        // Clean up
        delete refName;        
        
        return validDA;
    }
    
    void processJunctionVicinity(BamReader& reader, int32_t refLength, int32_t meanQueryLength, int32_t maxQueryLength, bool strandSpecific) {
        
        int32_t refId = intron->refId;
        Strand strand = intron->strand;

        uint32_t nbLeftFlankingAlignments = 0, nbRightFlankingAlignments = 0;

        int32_t regionStart = leftFlankStart - maxQueryLength - 1; regionStart < 0 ? 0 : regionStart;
        int32_t regionEnd = rightFlankEnd + maxQueryLength + 1; regionEnd >= refLength ? refLength - 1 : regionEnd;
        BamRegion region(refId, regionStart, refId, regionEnd);
        
        BamAlignment ba;
        
        // Focus only on the (expanded... to be safe...) region of interest
        if (!reader.SetRegion(region)) {
            BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                    "Could not set region")));
        }
        
        while(reader.GetNextAlignment(ba)) {
            
            // Look for left flanking alignments
            if (    intron->start > ba.Position && 
                    leftFlankStart <= ba.Position + ba.AlignedBases.size() &&
                    (!strandSpecific || strand == strandFromBool(ba.IsReverseStrand()))) {
                nbLeftFlankingAlignments++;
            }
            
            // Look for right flanking alignments
            if (    rightFlankEnd >= ba.Position && 
                    intron->end < ba.Position &&
                    (!strandSpecific || strand == strandFromBool(ba.IsReverseStrand()))) {
                nbRightFlankingAlignments++;
            }
        }
        
        this->setFlankingAlignmentCounts(nbLeftFlankingAlignments, nbRightFlankingAlignments);
    }
    
    /**
     * Call this method to recalculate all junction metrics based on the current location
     * and alignment information present in this junction
     * 
     * @param readLength
     */
    void calcAllRemainingMetrics(double readLength) {
       
        calcAnchorStats();      // Metrics 5 and 7
        calcEntropy();          // Metric 6
        calcAlignmentStats();   // Metrics 8 and 9
        calcMaxMMES();          // Metric 12
        calcMultipleMappingScore(); // Metric 18
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
     * Calculates the entropy score for this junction.  Higher entropy is generally
     * more indicative of a genuine junction than a lower score.
     *
     * This version of the function collects BamAlignment start positions and passes
     * them to an overloaded version of this function.
     * 
     * @return The entropy of this junction
     */
    double calcEntropy() {
        
        vector<int32_t> junctionPositions;
        
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            junctionPositions.push_back(ba.Position);
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
    double calcEntropy(vector<int32_t>& junctionPositions) {
        
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
    
    /**
     * Metric 13 and 14: Calculates the 5' and 3' hamming distances from a genomic
     * region represented by this junction
     * @param junctionSeq The DNA sequence representing this junction on the genome
     */
    void calcHammingScores(string& junctionSeq) {
        
        const uint8_t REGION_LENGTH = 10;
        string intron5p;
        string intron3p;
        string anchor5p;
        string anchor3p;
        
        int32_t intronStartOffset = intron->start - leftFlankStart;
        int32_t intronEndOffset = intron->end - leftFlankStart;
        
        // This test might be revisiting in the future.  It's possible that we could
        // just get a larger flanking region to calculate the many cases
        if (intron->size() < REGION_LENGTH 
                || intronStartOffset - REGION_LENGTH < 0 
                || intronEndOffset + REGION_LENGTH >= rightFlankEnd - leftFlankStart) {
            hammingDistance5p = -1;
            hammingDistance3p = -1;
        }
        else {

            Strand s = intron->strand != UNKNOWN ? intron->strand : predictedStrand;                
            
            switch(s) {
                case POSITIVE:
                    intron5p = junctionSeq.substr(intronStartOffset, REGION_LENGTH);
                    intron3p = junctionSeq.substr(intronEndOffset - REGION_LENGTH, REGION_LENGTH);
                    anchor5p = junctionSeq.substr(intronStartOffset - REGION_LENGTH, REGION_LENGTH);
                    anchor3p = junctionSeq.substr(intronEndOffset, REGION_LENGTH);
                    break;
                case NEGATIVE:
                    intron5p = SeqUtils::reverseComplement(junctionSeq.substr(intronEndOffset - REGION_LENGTH, REGION_LENGTH));
                    intron3p = SeqUtils::reverseComplement(junctionSeq.substr(intronStartOffset, REGION_LENGTH));
                    anchor5p = SeqUtils::reverseComplement(junctionSeq.substr(intronEndOffset, REGION_LENGTH));
                    anchor3p = SeqUtils::reverseComplement(junctionSeq.substr(intronStartOffset - REGION_LENGTH, REGION_LENGTH));
                    break;
                default:
                    hammingDistance5p = -1;
                    hammingDistance3p = -1;
                    return;
            }
            
            hammingDistance5p = SeqUtils::hammingDistance(anchor5p, intron3p);
            hammingDistance3p = SeqUtils::hammingDistance(anchor3p, intron5p);  
        }
    }
    
    /**
     * Calculates metric 12.  MaxMMES.
     */
    void calcMaxMMES() {
        
        uint16_t maxmmes = 0;
        
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            
            uint16_t leftMM = calcMinimalMatchInCigarDataSubset(ba, leftFlankStart, intron->start);
            uint16_t rightMM = calcMinimalMatchInCigarDataSubset(ba, intron->end, rightFlankEnd);
            
            uint16_t mmes = min(leftMM, rightMM);

            maxmmes = max(maxmmes, mmes);
        }
        
        this->maxMMES = maxmmes;
    }
    
    /**
     * Calculates metric 18.  Multiple mapping score
     */
    void calcMultipleMappingScore() {
        
        size_t N = this->getNbJunctionAlignments();
        
        uint32_t M = 0;
        BOOST_FOREACH(BamAlignment ba, junctionAlignments) {
            
            M += 1;  // Number of multiple splitting patterns
        }
        
        this->multipleMappingScore = (double)N / (double)M;
    }
    
    uint16_t calcMinimalMatchInCigarDataSubset(BamAlignment& ba, int32_t start, int32_t end) {
        
        if (start > ba.Position + ba.AlignedBases.size() || end < ba.Position)
            BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                    "Found an alignment that does not have a presence in the requested region")));
        
        int32_t pos = ba.Position;
        uint16_t mismatches = 0;
        int32_t length = 0;
        bool inRegion = false;
        
        BOOST_FOREACH(CigarOp op, ba.CigarData) {
           
            if (pos > end) {
                break;
            }
            
            inRegion = pos >= start;                
            
            if (BamUtils::opFollowsReference(op.Type)) {
                
                if (inRegion)
                    length += op.Length;

                pos += op.Length;
            }
            
            if (inRegion && op.Type == 'X') {
                mismatches++;
            }
        } 
        
        return length - mismatches;
    }
    
    double calcCoverage(int32_t a, int32_t b, const vector<uint32_t>& coverageLevels) {
        
        double multiplier = 1.0 / (b - a);
        uint32_t readCount = 0;
        
        for (int32_t i = a; i <= b; i++) {
            
            //int32_t pos = intron->strand == NEGATIVE ? i*2+1 : i*2;
            readCount += coverageLevels[i];
        }
        return multiplier * (double)readCount;
    }
    
    void calcCoverage(int32_t meanReadLength, const vector<uint32_t>& coverageLevels) {
        
        int32_t donorStart = intron->start - 2 * meanReadLength; donorStart = donorStart < 0 ? 0 : donorStart;
        int32_t donorMid = intron->start - meanReadLength;  // This one should be fine
        int32_t donorEnd = intron->start;                   // This is definitely ok, no capping required
        
        int32_t acceptorStart = intron->end - meanReadLength;
        int32_t acceptorMid = intron->end;
        int32_t acceptorEnd = intron->end + meanReadLength;
        
        
        if (donorMid - donorStart != donorEnd - donorMid ||
            acceptorMid - acceptorStart != acceptorEnd - acceptorMid) {
            
            coverage = -1.0;
        }
        else {
            
            double donorCoverage = 
                    calcCoverage(donorStart, donorMid-1, coverageLevels) -
                    calcCoverage(donorMid, donorEnd, coverageLevels);

            double acceptorCoverage = 
                    calcCoverage(acceptorMid, acceptorEnd, coverageLevels) -
                    calcCoverage(acceptorStart, acceptorMid-1, coverageLevels);

            coverage = donorCoverage + acceptorCoverage;
        }
    }
    
    
    
    // **** Core property getters ****
    
    shared_ptr<Intron> getIntron() const {
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
    bool hasCanonicalSpliceSites() const {
        return this->canonicalSpliceSites;
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

    /**
     * Metric 12: The maximum of the minimal mismatches exon sequences
     * @return 
     */
    uint32_t getMaxMMES() const {
        return maxMMES;
    }
    
    /**
     * Metric 13: Hamming distance between the 
     * @return 
     */
    int16_t getHammingDistance3p() const {
        return hammingDistance3p;
    }

    /**
     * Metric 14:
     * @return 
     */
    int16_t getHammingDistance5p() const {
        return hammingDistance5p;
    }

    /**
     * Metric 15:
     * @return 
     */
    double getCoverage() const {
        return coverage;
    }
    
    /**
     * Metric 16:
     * @return 
     */
    bool isUniqueJunction() const {
        return uniqueJunction;
    }

    /**
     * Metric 17:
     * @return 
     */
    bool isPrimaryJunction() const {
        return primaryJunction;
    }
    
    double getMultipleMappingScore() const {
        return multipleMappingScore;
    }

    Strand getPredictedStrand() const {
        return predictedStrand;
    }


    void setDa1(string da1) {
        this->da1 = da1;
    }

    void setDa2(string da2) {
        this->da2 = da2;
    }
    
    void setCoverage(double coverage) {
        this->coverage = coverage;
    }

    void setDiffAnchor(int32_t diffAnchor) {
        this->diffAnchor = diffAnchor;
    }

    void setEntropy(double entropy) {
        this->entropy = entropy;
    }

    void setHammingDistance3p(int16_t hammingDistance3p) {
        this->hammingDistance3p = hammingDistance3p;
    }

    void setHammingDistance5p(int16_t hammingDistance5p) {
        this->hammingDistance5p = hammingDistance5p;
    }

    void setMaxMMES(uint32_t maxMMES) {
        this->maxMMES = maxMMES;
    }

    void setMaxMinAnchor(int32_t maxMinAnchor) {
        this->maxMinAnchor = maxMinAnchor;
    }

    void setNbDistinctAlignments(uint32_t nbDistinctAlignments) {
        this->nbDistinctAlignments = nbDistinctAlignments;
    }

    void setNbDistinctAnchors(uint32_t nbDistinctAnchors) {
        this->nbDistinctAnchors = nbDistinctAnchors;
    }

    void setNbDownstreamFlankingAlignments(uint32_t nbDownstreamFlankingAlignments) {
        this->nbDownstreamFlankingAlignments = nbDownstreamFlankingAlignments;
    }

    void setNbReliableAlignments(uint32_t nbReliableAlignments) {
        this->nbReliableAlignments = nbReliableAlignments;
    }

    void setNbUpstreamFlankingAlignments(uint32_t nbUpstreamFlankingAlignments) {
        this->nbUpstreamFlankingAlignments = nbUpstreamFlankingAlignments;
    }
    
    void setMultipleMappingScore(double multipleMappingScore) {
        this->multipleMappingScore = multipleMappingScore;
    }

    void setPredictedStrand(Strand predictedStrand) {
        this->predictedStrand = predictedStrand;
    }


    
    // **** Output methods ****
    
    /**
     * Complete human readable description of this junction
     * @param strm
     */
    void outputDescription(std::ostream &strm) {
        outputDescription(strm, "\n");
    }
    /**
     * Complete human readable description of this junction
     * @param strm
     */
    void outputDescription(std::ostream &strm, string delimiter) {
        strm << "Location: ";
        
        if (intron != NULL) {
            intron->outputDescription(strm);
        }
        else {
            strm << "No location set";
        }
        
        strm << delimiter
             << "Flank limits: (" << leftFlankStart << ", " << rightFlankEnd << ")" << delimiter
             << "Junction Predictions:" << delimiter
             << "1:  Strand: " << strandToString(predictedStrand) << delimiter
             << "Junction Metrics:" << delimiter
             << "1:  # Junction Alignments: " << getNbJunctionAlignments() << delimiter
             << "2:  Canonical?: " << boolalpha << cssToString(canonicalSpliceSites) << "; Sequences: (" << da1 << " " << da2 << ")" << delimiter
             << "3:  Intron Size: " << getIntronSize() << delimiter
             << "4:  MaxMinAnchor: " << maxMinAnchor << delimiter
             << "5:  DiffAnchor: " << diffAnchor << delimiter
             << "6:  Entropy: " << entropy << delimiter
             << "7:  # Distinct Anchors: " << nbDistinctAnchors << delimiter
             << "8:  # Distinct Alignments: " << nbDistinctAlignments << delimiter
             << "9:  # Reliable (MP >=" << MAP_QUALITY_THRESHOLD << ") Alignments: " << nbReliableAlignments << delimiter
             << "10: # Upstream Non-Spliced Alignments: " << nbUpstreamFlankingAlignments << delimiter
             << "11: # Downstream Non-Spliced Alignments: " << nbDownstreamFlankingAlignments << delimiter
             << "12: MaxMMES: " << maxMMES << delimiter
             << "13: Hamming Distance 5': " << hammingDistance5p << delimiter
             << "14: Hamming Distance 3': " << hammingDistance3p << delimiter
             << "15: Coverage: " << coverage << delimiter
             << "16: Unique Junction: " << boolalpha << uniqueJunction << delimiter
             << "17: Primary Junction: " << boolalpha << primaryJunction << delimiter
             << "18: Multiple mapping score: " << multipleMappingScore << delimiter;
    }
    
    /**
     * Complete human readable description of this junction
     * @param strm
     */
    void outputGFF(std::ostream &strm, uint32_t id, RefVector& refs) {
        
        // Use intron strand if known, otherwise use the predicted strand,
        // if predicted strand is also unknown then use "." to indicated unstranded
        const char strand = (intron->strand == UNKNOWN ?
                                predictedStrand == UNKNOWN ?
                                    '.' :
                                    strandToChar(predictedStrand) :
                                strandToChar(intron->strand));
        
        string juncId = string("junc_") + lexical_cast<string>(id);
        string refName = refs[intron->refId].RefName;
        
        // Output junction parent
        strm << refName << "\t"
             << "portculis" << "\t"     // source
             << "junction" << "\t"      // type (may change later)
             << leftFlankStart << "\t"  // start
             << rightFlankEnd << "\t"   // end
             << "0.0" << "\t"           // No score for the moment
             << strand << "\t"          // strand
             << "." << "\t"             // Just put "." for the phase
             << "ID=" << juncId << ";"
             << "Note=";
        outputDescription(strm, "%0A");
        strm << endl;

        // Output left exonic region
        strm << refName << "\t"
             << "portculis" << "\t"
             << "partial_exon" << "\t"
             << leftFlankStart << "\t"
             << (intron->start - 1) << "\t"
             << "0.0" << "\t"
             << strand << "\t"
             << "." << "\t"
             << "ID=" << juncId << "_left" << ";"
             << "Parent=" << juncId << endl;
                
        // Output right exonic region
        strm << refName << "\t"
             << "portculis" << "\t"
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
    void outputBED(std::ostream &strm, uint32_t id, RefVector& refs) {
        
        // Use intron strand if known, otherwise use the predicted strand,
        // if predicted strand is also unknown then use "." to indicated unstranded
        const char strand = (intron->strand == UNKNOWN ?
                                predictedStrand == UNKNOWN ?
                                    '.' :
                                    strandToChar(predictedStrand) :
                                strandToChar(intron->strand));
        
        string juncId = string("junc_") + lexical_cast<string>(id);
        string refName = refs[intron->refId].RefName;
        
        int32_t sz1 = intron->start - leftFlankStart;
        int32_t sz2 = rightFlankEnd - intron->end;
        string blockSizes = lexical_cast<string>(sz1) + "," + lexical_cast<string>(sz2);
        string blockStarts = lexical_cast<string>(0) + "," + lexical_cast<string>(intron->end+1 - leftFlankStart);
        
        // Output junction parent
        strm << refName << "\t"         // chrom
             << leftFlankStart << "\t"  // chromstart
             << rightFlankEnd << "\t"   // chromend
             << juncId << "\t"          // name
             << "0.0" << "\t"           // No score for the moment
             << strand << "\t"          // strand
             << intron->start << "\t"   // thickstart
             << intron->end << "\t"     // thickend
             << "255,0,0" << "\t"       // Just use red for the moment
             << "2" << "\t"             // 2 blocks: Left and right block
             << blockSizes << "\t"
             << blockStarts << endl;        
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
                    << j.da1 << "\t"
                    << j.da2 << "\t"
                    << strandToChar(j.predictedStrand) << "\t"
                    << j.getNbJunctionAlignments() << "\t"
                    << cssToChar(j.canonicalSpliceSites) << "\t"
                    << j.getIntronSize() << "\t"
                    << j.maxMinAnchor << "\t"
                    << j.diffAnchor << "\t"
                    << j.entropy << "\t"
                    << j.nbDistinctAnchors << "\t"
                    << j.nbDistinctAlignments << "\t"
                    << j.nbReliableAlignments << "\t"
                    << j.nbUpstreamFlankingAlignments << "\t"
                    << j.nbDownstreamFlankingAlignments << "\t"
                    << j.maxMMES << "\t"
                    << j.hammingDistance5p << "\t"
                    << j.hammingDistance3p << "\t"
                    << j.coverage << "\t"
                    << j.uniqueJunction << "\t"
                    << j.primaryJunction << "\t"
                    << j.multipleMappingScore;
    }
    
        
    // **** Static methods ****

    /**
     * Header for table output
     * @return 
     */
    static string junctionOutputHeader() {
        return string(Intron::locationOutputHeader()) + "\tleft\tright\tss1\tss2\t" + 
                boost::algorithm::join(PREDICTION_NAMES, "\t") + "\t" +
                boost::algorithm::join(METRIC_NAMES, "\t");
    }

    static shared_ptr<Junction> parse(const string& line) {
        
        vector<string> parts; // #2: Search for tokens
        boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

        if (parts.size() != 24) {
            BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
                "Could not parse line due to incorrect number of columns. Expected 24 columns: ") + line));
        }

        // Create intron
        shared_ptr<Intron> i(new Intron(
            lexical_cast<int32_t>(parts[1]),
            lexical_cast<int32_t>(parts[2]),
            lexical_cast<int32_t>(parts[3]),
            strandFromChar(parts[4][0])
        ));

        // Create basic junction
        shared_ptr<Junction> j(new Junction(
            i,
            lexical_cast<int32_t>(parts[5]),
            lexical_cast<int32_t>(parts[6])
        ));

        // Splice site strings
        j->setDa1(parts[7]);
        j->setDa2(parts[8]);

        // Set predictions to junction
        j->setPredictedStrand(strandFromChar(parts[9][0]));
                
        // Add metrics to junction
        j->setDonorAndAcceptorMotif(cssFromChar(parts[11][0]));
        j->setMaxMinAnchor(lexical_cast<int32_t>(parts[13]));
        j->setDiffAnchor(lexical_cast<int32_t>(parts[14]));
        j->setEntropy(lexical_cast<double>(parts[15]));
        j->setNbDistinctAnchors(lexical_cast<uint32_t>(parts[16]));
        j->setNbDistinctAlignments(lexical_cast<uint32_t>(parts[17]));
        j->setNbReliableAlignments(lexical_cast<uint32_t>(parts[18]));
        j->setNbUpstreamFlankingAlignments(lexical_cast<uint32_t>(parts[19]));
        j->setNbDownstreamFlankingAlignments(lexical_cast<uint32_t>(parts[20]));
        j->setMaxMMES(lexical_cast<uint32_t>(parts[21]));
        j->setHammingDistance5p(lexical_cast<uint32_t>(parts[22]));
        j->setHammingDistance3p(lexical_cast<uint32_t>(parts[23]));
        j->setCoverage(lexical_cast<double>(parts[24]));
        j->setUniqueJunction(lexical_cast<bool>(parts[25]));
        j->setPrimaryJunction(lexical_cast<bool>(parts[26]));
        j->setMultipleMappingScore(lexical_cast<double>(parts[27]));

        
        return j;
    }
};

}
