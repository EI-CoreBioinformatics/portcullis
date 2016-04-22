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

#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
using std::ostream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::string;
using std::size_t;
using std::vector;
using std::shared_ptr;

typedef std::unordered_map<size_t, uint16_t> SplicedAlignmentMap;

#include <boost/exception/all.hpp>

#include "bam/bam_master.hpp"
#include "bam/bam_alignment.hpp"
#include "bam/bam_reader.hpp"
#include "bam/genome_mapper.hpp"
using namespace portcullis::bam;

#include "intron.hpp"
#include "seq_utils.hpp"
using portcullis::Intron;


namespace portcullis {    

const uint16_t MAP_QUALITY_THRESHOLD = 30;

typedef boost::error_info<struct JunctionError,string> JunctionErrorInfo;
struct JunctionException: virtual boost::exception, virtual std::exception { };

const vector<string> JO_NAMES = {
        "JO01",
        "JO02",
        "JO03",
        "JO04",
        "JO05",
        "JO06",
        "JO07",
        "JO08",
        "JO09",
        "JO10",
        "JO11",
        "JO12",
        "JO13",
        "JO14",
        "JO15",
        "JO16",
        "JO17",
        "JO18",
        "JO19",
        "JO20"
};

const vector<string> METRIC_NAMES = {
        "M1-canonical_ss",
        "M2-nb_reads",
        "M3-nb_dist_aln",
        "M4-nb_rel_aln",
        "M5-intron_size",
        "M6-left_anc_size",
        "M7-right_anc_size",
        "M8-max_min_anc",
        "M9-dif_anc",
        "M10-dist_anc",
        "M11-entropy",
        "M12-maxmmes",
        "M13-hamming5p",
        "M14-hamming3p",
        "M15-coverage",
        "M16-uniq_junc",
        "M17-primary_junc",
        "M18-mm_score",
        "M19-mean_mismatches",
        "M20-nb_usrs",
        "M21-nb_msrs",
        "M22-rel2raw",
        "M23-nb_up_juncs",
        "M24-nb_down_juncs",
        "M25-up_aln",
        "M26-down_aln",
        "M27-dist_2_up_junc",
        "M28-dist_2_down_junc",
        "M29-dist_nearest_junc"
    };

const vector<string> STRAND_NAMES = {
        "read-strand",
        "ss-strand",
        "consensus-strand"
    };

const string CANONICAL_SEQ = "GTAG";
const string SEMI_CANONICAL_SEQ1 = "ATAC";
const string SEMI_CANONICAL_SEQ2 = "GCAG";

const string CANONICAL_SEQ_RC = SeqUtils::reverseComplement(CANONICAL_SEQ);
const string SEMI_CANONICAL_SEQ1_RC = SeqUtils::reverseComplement(SEMI_CANONICAL_SEQ1);
const string SEMI_CANONICAL_SEQ2_RC = SeqUtils::reverseComplement(SEMI_CANONICAL_SEQ2);

// Represents both upstream and downstream coverage levels.  The half way point represents
// the junction.  Use an even number to ensure upstream and downstream lengths are
// equal and to avoid any issues.
const uint32_t TRIMMED_COVERAGE_LENGTH = 50;

enum CanonicalSS {
    CANONICAL,
    SEMI_CANONICAL,
    NO,
    ALL
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



struct AlignmentInfo {
    BamAlignmentPtr ba;
    size_t nameCode;
    uint32_t totalUpstreamMatches; // Total number of upstream matches in this junction window
    uint32_t totalDownstreamMatches; // Total number of downstream matches in this junction window
    uint32_t totalUpstreamMismatches;
    uint32_t totalDownstreamMismatches;
    uint32_t upstreamMatches;   // Distance to first upstream mismatch
    uint32_t downstreamMatches; // Distance to first downstream mismatch
    uint32_t minMatch;  // Distance to first mismatch (minimum of either upstream or downstream)
    uint32_t maxMatch;  // Distance to first mismatch (maximum of either upstream or downstream)
    uint32_t nbMismatches; // Total number of mismatches in this junction window
    uint32_t mmes; // Minimal Match on Either Side of exon junction    
    
    AlignmentInfo(BamAlignmentPtr _ba) {
        
        // Copy alignment
        ba = _ba;
        
        // Calculate a hash of the alignment name
        nameCode = std::hash<std::string>()(ba->deriveName());
        
        totalUpstreamMatches = 0;
        totalDownstreamMatches = 0;
        totalUpstreamMismatches = 0;
        totalDownstreamMismatches = 0;
        upstreamMatches = 0;
        downstreamMatches = 0;
        minMatch = 0;
        maxMatch = 0;
        nbMismatches = 0;
        mmes = 0;
    }
    
    ~AlignmentInfo() {        
        ba.reset();       
    }
    
    void calcMatchStats(const Intron& i, const uint32_t leftStart, const uint32_t rightEnd, const string& ancLeft, const string& ancRight);
 
    uint32_t getNbMatchesFromStart(const string& query, const string& anchor);
    uint32_t getNbMatchesFromEnd(const string& query, const string& anchor);
};

typedef shared_ptr<AlignmentInfo> AlignmentInfoPtr;

class Junction {
    
    
private:
    
    // **** Properties that describe where the junction is ****
    shared_ptr<Intron> intron;
    vector<shared_ptr<AlignmentInfo>> alignments;
    vector<size_t> alignmentCodes;
    vector<uint32_t> trimmedCoverage;
    vector<double> trimmedLogDevCov;
    uint32_t meanQueryLength;
    bool suspicious;
    bool pfp;
    
    // **** Junction metrics ****
    CanonicalSS canonicalSpliceSites;           // Metric 1
    uint32_t nbJunctionAlignments;              // Metric 2
    uint32_t nbDistinctAlignments;              // Metric 3
    uint32_t nbReliableAlignments;              // Metric 4
                                                // Metric 5 (intron size) calculated via location properties
    uint32_t leftAncSize;                       // Metric 6
    uint32_t rightAncSize;                      // Metric 7
    int32_t  maxMinAnchor;                      // Metric 8
    int32_t  diffAnchor;                        // Metric 9
    uint32_t nbDistinctAnchors;                 // Metric 10
    double   entropy;                           // Metric 11
    uint32_t maxMMES;                           // Metric 12
    int16_t  hammingDistance5p;                 // Metric 13
    int16_t  hammingDistance3p;                 // Metric 14
    double   coverage;                          // Metric 15
    bool     uniqueJunction;                    // Metric 16
    bool     primaryJunction;                   // Metric 17    
    double   multipleMappingScore;              // Metric 18
    double   meanMismatches;                    // Metric 19
    //uint32_t nbUniquelySplicedReads;          // Metric 20 (Use getter)
    uint32_t nbMultipleSplicedReads;            // Metric 21
    //double   reliableVsRawReadRatio;            // Metric 22 (Use getter)
    uint16_t nbUpstreamJunctions;               // Metric 23
    uint16_t nbDownstreamJunctions;             // Metric 24
    uint32_t nbUpstreamFlankingAlignments;      // Metric 25
    uint32_t nbDownstreamFlankingAlignments;    // Metric 26
    int32_t distanceToNextUpstreamJunction;     // Metric 27
    int32_t distanceToNextDownstreamJunction;   // Metric 28
    int32_t distanceToNearestJunction;          // Metric 29
    vector<uint32_t> junctionOverhangs;         // Metric 30-49
    
    // **** Predictions ****
    
    Strand readStrand;      // Strand derived from alignments
    Strand ssStrand;        // Strand derived from splice sites    
    Strand consensusStrand; // If readStrand and ssStrand agree then strand is the same, otherwise UNKNOWN
    
    bool genuine;           // Used as a hidden variable for use with cross validating a trained model instance.
    
    // **** Additional properties ****
    
    int32_t leftAncStart;
    int32_t rightAncEnd;
    string da1, da2;                    // These store the dinucleotides found at the predicted donor / acceptor sites in the intron
    
    
protected:
    
    /**
     * Tests whether the two strings could represent valid donor and acceptor sites
     * for this junction
     * @param seq1
     * @param seq2
     * @return 
     */
    CanonicalSS hasCanonicalSpliceSites(const string& seq1, const string& seq2);
    
    Strand predictedStrandFromSpliceSites(const string& seq1, const string& seq2);
    
    
public:
    
    // **** Constructors ****
    
    Junction(shared_ptr<Intron> _location, int32_t _leftAncStart, int32_t _rightAncEnd);
    
    /**
     * Copy constructor
     * @param j The other junction to deep copy into this
     */
    Junction(const Junction& j) : Junction(j, true) {}
        
    /**
     * Copy constructor, with option to copy bam alignments associated with the junction
     * @param j The other junction to deep copy into this
     * @param withAlignments Whether to copy over the alignments or not
     */
    Junction(const Junction& j, bool withAlignments);
    
    // **** Destructor ****
    virtual ~Junction();
    
    void clearAlignments();
    
    const BamAlignment& getFirstAlignment() const {
        return *(alignments[0]->ba);
    }
    
   
    const Intron& getLocation() const {
        return *intron;
    }

    void setLocation(shared_ptr<Intron> location) {
        this->intron = location;
    }

    uint32_t getMeanQueryLength() const {
        return meanQueryLength;
    }

    void setMeanQueryLength(uint32_t meanQueryLength) {
        this->meanQueryLength = meanQueryLength;
    }

    
    
    void addJunctionAlignment(const BamAlignment& al);
    
    void setNbJunctionAlignments(uint32_t nbJunctionAlignments) {
        this->nbJunctionAlignments = nbJunctionAlignments;
    }
    
    void setDonorAndAcceptorMotif(CanonicalSS canonicalSS) {
        this->canonicalSpliceSites = canonicalSS;
    }
    
    CanonicalSS setDonorAndAcceptorMotif(string seq1, string seq2);
    
    void setFlankingAlignmentCounts(uint32_t nbUpstream, uint32_t nbDownstream);
    
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
     * Extends the anchor regions of the junction.  Also updates any relevant 
     * metrics.
     * @param otherStart The alternative start position of the left anchor
     * @param otherEnd The alternative end position of the right anchor
     */
    void extendAnchors(int32_t otherStart, int32_t otherEnd);
        
    void determineStrandFromReads();
    
    void processJunctionWindow(const GenomeMapper& genomeMapper);
    
    void processJunctionVicinity(BamReader& reader, int32_t refLength, int32_t meanQueryLength, int32_t maxQueryLength);
    
    
    void calcMetrics();
    
    void updateAlignmentInfo();
    
    /**
     * Metric 5 and 7: Diff Anchor and # Distinct Anchors
     * @return 
     */
    void calcAnchorStats();
    
    
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
     * @return The entropy of this junction
     */
    double calcEntropy();
    
    /**
     * Provides a convienient way of get the shannon entropy score outside the
     * junction / junction system framework
     * @param offsets
     * @return 
     */
    double calcEntropy(const vector<int32_t> offsets);
    
    /**
     * Metrics: # Distinct Alignments, # Unique/Reliable Alignments, #mismatches
     * @return 
     */
    void calcAlignmentStats();
    
    /**
     * Metric 13 and 14: Calculates the 5' and 3' hamming distances from a genomic
     * region represented by this junction
     */
    void calcHammingScores( const string& leftAnchor, const string& leftIntron, 
                            const string& rightIntron, const string& rightAnchor);
    
    /**
     * Calculates MaxMMES, mismatches, and junction overhangs.
     * Requires alignment information to be populated
     */
    void calcMismatchStats();
        
    /**
     * Calculates metric 18.  Multiple mapping score
     */
    void calcMultipleMappingScore(SplicedAlignmentMap& map);
    

    void calcTrimmedCoverageVector();
    
    double calcCoverage(int32_t a, int32_t b, const vector<uint32_t>& coverageLevels);
    
    double calcCoverage(const vector<uint32_t>& coverageLevels);
    
    
    
    // **** Core property getters ****
    
    shared_ptr<Intron> getIntron() const {
        return intron;
    }

    int32_t getLeftAncStart() const {
        return leftAncStart;
    }

    int32_t getRightAncEnd() const {
        return rightAncEnd;
    }
    
    size_t size() const {
        return rightAncEnd - leftAncStart + 1;
    }

    
    size_t getNbJunctionAlignmentFromVector() const {
        return this->alignments.size();
    }
    
    bool isSuspicious() const {
        return this->suspicious;
    }
    
    void setSuspicious(bool suspicious) {
        this->suspicious = suspicious;
    }
    
    bool isPotentialFalsePositive() const {
        return this->pfp;
    }
    
    void setPotentialFalsePositive(bool pfp) {
        this->pfp = pfp;
    }
    
    // **** Metric getters ****
    
    /**
     * The number of alignments directly supporting this junction
     * @return 
     */
    size_t getNbJunctionAlignments() const {
        return this->nbJunctionAlignments;
    }
    
    /**
     * Whether or not there is a donor and acceptor motif at the two base
     * pairs at the start and end of the junction / intron
     * @return 
     */
    bool hasCanonicalSpliceSites() const {
        return this->canonicalSpliceSites;
    }
    
    CanonicalSS getSpliceSiteType() const {
        return this->canonicalSpliceSites;
    }
    
    /**
     * The intron size
     * @return 
     */
    int32_t getIntronSize() const {
        return intron != nullptr ? intron->size() : 0;
    }
    
    int32_t getLeftAnchorSize() const {
        return leftAncSize;
    }
    
    int32_t getRightAnchorSize() const {
        return rightAncSize;                
    }
    
    /**
     * The maximum anchor distance from the shortest side of each supporting 
     * alignment
     * @return
     */
    int32_t getMaxMinAnchor() const {
        return this->maxMinAnchor;
    }
    
    /**
     * Diff Anchor
     * @return 
     */
    int32_t getDiffAnchor() const {        
        return this->diffAnchor;
    }
    
    /**
     * Entropy (definition from "Graveley et al, The developmental 
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
    
    
    uint32_t getNbDistinctAnchors() const {
        return nbDistinctAnchors;
    }
    
    uint32_t getNbDistinctAlignments() const {
        return nbDistinctAlignments;
    }

    uint32_t getNbReliableAlignments() const {
        return nbReliableAlignments;
    }

    /**
     * The number of upstream non-spliced supporting reads
     * @return 
     */
    uint32_t getNbUpstreamFlankingAlignments() const {
        return nbUpstreamFlankingAlignments;
    }
    
    /**
     * The number of downstream non-spliced supporting reads
     * @return 
     */
    uint32_t getNbDownstreamFlankingAlignments() const {
        return nbDownstreamFlankingAlignments;
    }

    /**
     * The maximum of the minimal mismatches exon sequences
     * @return 
     */
    uint32_t getMaxMMES() const {
        return maxMMES;
    }
    
    /**
     * Hamming distance between the 
     * @return 
     */
    int16_t getHammingDistance3p() const {
        return hammingDistance3p;
    }

    int16_t getHammingDistance5p() const {
        return hammingDistance5p;
    }

    /**
     * Metric 15: The coverage score around the junction.  Uses unspliced reads
     * around both donor and acceptor sites.  If this is geniune junction you expect
     * unspliced reads to drop off closer to the donor and acceptor sites you get
     * (relative to read length).  We compare a window one read away from the d/a
     * site and one read up to the d/a site.  We then substract the mean coverage from
     * the more distant window to the close window.  We then add these score from
     * both donor and acceptor sites together.  For genuine junctions this should be
     * the score should be relatively large.  See TrueSight paper for more info.
     * @return 
     */
    double getCoverage() const {
        return coverage;
    }
    
    bool isUniqueJunction() const {
        return uniqueJunction;
    }

    bool isPrimaryJunction() const {
        return primaryJunction;
    }
    
    /**
     * Metric 18: Multiple mapping score (reflects mapping ambiguity.  Small score
     * implies reads in this junction align to other splice junctions across the
     * genome).  From TrueSight paper.
     * @return 
     */
    double getMultipleMappingScore() const {
        return multipleMappingScore;
    }
    
    /**
     * Metric 19: Number of mismatches.  The total number of mismatches from all
     * splice aligned reads in this junction.  From TrueSight paper.
     * @return 
     */
    double getMeanMismatches() const {
        return meanMismatches;
    }
    
    /**
     * The number of spliced reads in this junction that also cover additional junctions
     * @return 
     */
    uint32_t getNbMultipleSplicedReads() const {
        return nbMultipleSplicedReads;
    }
    
    uint32_t getNbUniquelySplicedReads() const {
        return nbJunctionAlignments - nbMultipleSplicedReads;
    }
    
    double getReliable2RawRatio() const {
        return (double)nbReliableAlignments / (double)nbJunctionAlignments;
    }
    
    uint16_t getNbDownstreamJunctions() const {
        return nbDownstreamJunctions;
    }

    uint16_t getNbUpstreamJunctions() const {
        return nbUpstreamJunctions;
    }
    
    uint32_t getDistanceToNearestJunction() const {
        return distanceToNearestJunction;
    }

    void setDistanceToNearestJunction(uint32_t distanceToNearestJunction) {
        this->distanceToNearestJunction = distanceToNearestJunction;
    }

    uint32_t getDistanceToNextDownstreamJunction() const {
        return distanceToNextDownstreamJunction;
    }

    void setDistanceToNextDownstreamJunction(uint32_t distanceToNextDownstreamJunction) {
        this->distanceToNextDownstreamJunction = distanceToNextDownstreamJunction;
    }

    uint32_t getDistanceToNextUpstreamJunction() const {
        return distanceToNextUpstreamJunction;
    }

    void setDistanceToNextUpstreamJunction(uint32_t distanceToNextUpstreamJunction) {
        this->distanceToNextUpstreamJunction = distanceToNextUpstreamJunction;
    }

    Strand getReadStrand() const {
        return readStrand;
    }

    Strand getSpliceSiteStrand() const {
        return ssStrand;
    }

    Strand getConsensusStrand() const {
        return consensusStrand;
    }

    void setDa1(string da1) {
        this->da1 = da1;
    }

    void setDa2(string da2) {
        this->da2 = da2;
    }
    
    void setLeftAncSize(uint16_t leftAncSize) {
        this->leftAncSize = leftAncSize;
    }

    void setNbDownstreamJunctions(uint16_t nbDownstreamJunctions) {
        this->nbDownstreamJunctions = nbDownstreamJunctions;
    }

    void setNbUpstreamJunctions(uint16_t nbUpstreamJunctions) {
        this->nbUpstreamJunctions = nbUpstreamJunctions;
    }

    void setRightAncSize(uint16_t rightAncSize) {
        this->rightAncSize = rightAncSize;
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

    void setMeanMismatches(double meanMismatches) {
        this->meanMismatches = meanMismatches;
    }
        
    void setNbMultipleSplicedReads(uint32_t nbMultipleSplicedReads) {
        this->nbMultipleSplicedReads = nbMultipleSplicedReads;
    }

    void setTrimmedLogDevCov(vector<double> trimmedLogDevCov) {
        this->trimmedLogDevCov = trimmedLogDevCov;
    }
    
    uint32_t getJunctionOverhangs(size_t index) const {
        return junctionOverhangs[index];
    }
    
    void setJunctionOverhangs(size_t index, uint32_t val) {
        junctionOverhangs[index] = val;
    }


    bool isGenuine() const {
        return genuine;
    }

    void setGenuine(bool genuine) {
        this->genuine = genuine;
    }
    
    string locationAsString() const {
        return this->intron->toString() + strandToChar(this->consensusStrand);
    }
    
    /**
     * Calculates a score for this intron size based on how this intron size fits
     * into an expected distribution specified by the length at the threhsold percentile
     * (threshold) provided by the user.  Introns of length < threshold have a score of 0. 
     * Introns with length > threshold have score: -ln(size - length_threshold)
     * @param Length of threshold Intron size
     * @return A score for this intron size given the threshold value
     */
    double calcIntronScore(const uint32_t threshold) const {
        return this->intron->size() <= threshold ? 0.0 : log(this->intron->size() - threshold);
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
    void outputDescription(std::ostream &strm, string delimiter);
    
    /**
     * Key value pair representation... useful for GFF attributes
     * @param strm
     * @param delimiter
     */
    void condensedOutputDescription(std::ostream &strm, string delimiter);
    
    
    /**
     * Complete human readable description of this intron (for augustus hints)
     * @param strm
     */
    void outputIntronGFF(std::ostream &strm, uint32_t id, const string& source);
    
    /**
     * Complete human readable description of this junction
     * @param strm
     */
    void outputJunctionGFF(std::ostream &strm, uint32_t id, const string& source);
    
    /**
     * Complete human readable description of this junction
     * @param strm
     */
    void outputBED(std::ostream &strm, const string& prefix, uint32_t id);
    
    
    /**
     * Represents this junction as a table row
     * @param seq1
     * @param seq2
     * @return 
     */
    friend ostream& operator<<(ostream &strm, Junction& j) {
        strm << *(j.intron) << "\t"
                    << j.leftAncStart << "\t"
                    << j.rightAncEnd << "\t"
                    << j.da1 << "\t"
                    << j.da2 << "\t"
                    << strandToChar(j.readStrand) << "\t"
                    << strandToChar(j.ssStrand) << "\t"
                    << strandToChar(j.consensusStrand) << "\t"
                    << cssToChar(j.canonicalSpliceSites) << "\t"
                    << j.getNbJunctionAlignments() << "\t"
                    << j.nbDistinctAlignments << "\t"
                    << j.nbReliableAlignments << "\t"
                    << j.getIntronSize() << "\t"
                    << j.getLeftAnchorSize() << "\t"
                    << j.getRightAnchorSize() << "\t"
                    << j.maxMinAnchor << "\t"
                    << j.diffAnchor << "\t"
                    << j.nbDistinctAnchors << "\t"
                    << j.entropy << "\t"
                    << j.maxMMES << "\t"
                    << j.hammingDistance5p << "\t"
                    << j.hammingDistance3p << "\t"
                    << j.coverage << "\t"
                    << j.uniqueJunction << "\t"
                    << j.primaryJunction << "\t"
                    << j.multipleMappingScore << "\t"
                    << j.meanMismatches << "\t"
                    << j.getNbUniquelySplicedReads() << "\t"
                    << j.nbMultipleSplicedReads << "\t"
                    << j.getReliable2RawRatio() << "\t"
                    << j.nbDownstreamJunctions << "\t"
                    << j.nbUpstreamJunctions << "\t"
                    << j.nbUpstreamFlankingAlignments << "\t"
                    << j.nbDownstreamFlankingAlignments << "\t"
                    << j.distanceToNextUpstreamJunction << "\t"
                    << j.distanceToNextDownstreamJunction << "\t"
                    << j.distanceToNearestJunction << "\t"
                    << j.meanQueryLength << "\t"
                    << j.suspicious << "\t"
                    << j.pfp;
        
        for(size_t i = 0; i < JO_NAMES.size(); i++) {            
            strm << "\t" << j.junctionOverhangs[i];
        }

        return strm;
    }
    
        
    // **** Static methods ****

    /**
     * Header for table output
     * @return 
     */
    static string junctionOutputHeader();

    static shared_ptr<Junction> parse(const string& line);
};


struct JunctionComparator {
    inline bool operator() (const shared_ptr<Junction>& j1, const shared_ptr<Junction>& j2) const {
        return IntronComparator()(*(j1->getIntron()), *(j2->getIntron()));        
    }
};

typedef shared_ptr<Junction> JunctionPtr;
typedef std::vector<shared_ptr<Junction> > JunctionList;

}
