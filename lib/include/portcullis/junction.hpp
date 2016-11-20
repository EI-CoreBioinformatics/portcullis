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
#include <map>
using std::ostream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::string;
using std::size_t;
using std::vector;
using std::map;
using std::shared_ptr;

typedef std::unordered_map<size_t, uint16_t> SplicedAlignmentMap;

#include <boost/exception/all.hpp>

#include "bam/bam_master.hpp"
#include "bam/bam_alignment.hpp"
#include "bam/bam_reader.hpp"
#include "bam/genome_mapper.hpp"
using namespace portcullis::bam;

#include "ml/markov_model.hpp"
using portcullis::ml::KmerMarkovModel;
using portcullis::ml::PosMarkovModel;

#include "intron.hpp"
#include "seq_utils.hpp"
using portcullis::Intron;


namespace portcullis {

/**
 *  Value of 30 is selected as it appears to be a good cutoff for distinguishing 
 * alignments likely to be uniquely mapped from those that have a high probability of
 * mapping to multiple locations.  This cutoff seems to work for most RNAseq mappers
 * that we are aware of.
 */ 
const uint16_t MAP_QUALITY_THRESHOLD = 30;

typedef boost::error_info<struct JunctionError, string> JunctionErrorInfo;

struct JunctionException : virtual boost::exception, virtual std::exception {
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

enum class CanonicalSS {
	CANONICAL,
	SEMI_CANONICAL,
	NO,
	ALL
};

inline CanonicalSS cssFromChar(char css) {
	switch (css) {
	case 'C':
		return CanonicalSS::CANONICAL;
	case 'S':
		return CanonicalSS::SEMI_CANONICAL;
	case 'N':
		return CanonicalSS::NO;
	}

	return CanonicalSS::NO;
}

inline char cssToChar(CanonicalSS css) {

	switch (css) {
	case CanonicalSS::CANONICAL:
		return 'C';
	case CanonicalSS::SEMI_CANONICAL:
		return 'S';
	case CanonicalSS::NO:
		return 'N';
	}

	return 'N';
}

inline string cssToString(CanonicalSS css) {

	switch (css) {
	case CanonicalSS::CANONICAL:
		return "Canonical";
	case CanonicalSS::SEMI_CANONICAL:
		return "Semi-canonical";
	case CanonicalSS::NO:
		return "No";
	}

	return string("No");
}

struct SplicingScores {
	double positionWeighting = 0.0;
	double splicingSignal = 0.0;
};

struct AlignmentInfo {
	BamAlignmentPtr ba;
	size_t nameCode;
	uint32_t totalUpstreamMatches; // Total number of upstream matches in this junction window
	uint32_t totalDownstreamMatches; // Total number of downstream matches in this junction window
	uint32_t totalUpstreamMismatches;
	uint32_t totalDownstreamMismatches;
	uint32_t upstreamMatches; // Distance to first upstream mismatch
	uint32_t downstreamMatches; // Distance to first downstream mismatch
	uint32_t minMatch; // Distance to first mismatch (minimum of either upstream or downstream)
	uint32_t maxMatch; // Distance to first mismatch (maximum of either upstream or downstream)
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


	// **** Junction metrics ****

	CanonicalSS canonicalSpliceSites; 

	// Alignment counts
	uint32_t nbAlRaw;
	uint32_t nbAlDistinct;
	uint32_t nbAlMultiplySpliced; // uniquely spliced = split - multiply spliced
	uint32_t nbAlUniquelyMapped; // multiply mapped = split - uniquely_mapped
	uint32_t nbAlBamProperlyPaired;
	uint32_t nbAlPortcullisProperlyPaired;
	uint32_t nbAlReliable; // reliable to raw ratio = reliable / split

	// RNAseq Derived Stats
	double entropy;
	double meanMismatches;
	double meanReadLength;
	uint32_t maxMinAnchor;
	uint32_t maxMMES;
	double intronScore;

	// Genome properties
	int16_t hammingDistance5p; // Metric 16
	int16_t hammingDistance3p; // Metric 17

	// Junction relationships
	bool uniqueJunction; // Metric 18
	bool primaryJunction; // Metric 19        
	uint16_t nbDownstreamJunctions; // Metric 20
	uint16_t nbUpstreamJunctions; // Metric 21
	int32_t distanceToNextDownstreamJunction; // Metric 22
	int32_t distanceToNextUpstreamJunction; // Metric 23
	int32_t distanceToNearestJunction; // Metric 24

	// Metrics derived from extra processing
	double multipleMappingScore; // Metric 25
	double coverage; // Metric 26
	uint32_t nbDownstreamFlankingAlignments; // Metric 27
	uint32_t nbUpstreamFlankingAlignments; // Metric 28

	// Junction anchor depth metrics
	vector<uint32_t> trimmedCoverage;
	vector<double> trimmedLogDevCov;
	bool suspicious;
	bool pfp;
	vector<uint32_t> junctionAnchorDepth;


	// **** Predictions ****

	Strand readStrand; // Strand derived from alignments
	Strand ssStrand; // Strand derived from splice sites    
	Strand consensusStrand; // If readStrand and ssStrand agree then strand is the same, otherwise UNKNOWN
	double score; // Only applied if the random forest makes a prediction


	// **** Additional properties ****

	int32_t leftAncStart;
	int32_t rightAncEnd;
	string da1, da2; // These store the dinucleotides found at the predicted donor / acceptor sites in the intron    
	uint32_t id; // Unique identifier for the junction
	bool genuine; // Used as a hidden variable for use with cross validating a trained model instance.

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
	Junction(const Junction& j) : Junction(j, true) {
	}

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



	// ***** Getters *****

	/**
	 * The strand according to 95% of alignments in this junction.  Only
	 * used if strand-specific mode is provided by the user.  Otherwise unknown
	 * strand.
	 * @return 
	 */
	Strand getReadStrand() const {
		return readStrand;
	}

	/**
	 * The strand according to the splice sites of this junction.  Only set if
	 * the splice sites are canonical or semi-canonical
	 * @return 
	 */
	Strand getSpliceSiteStrand() const {
		return ssStrand;
	}

	/**
	 * The strand according to a consensus or read strand and splice site strand.
	 * If there is disagreement then this is set to unknown strand.
	 * @return 
	 */
	Strand getConsensusStrand() const {
		return consensusStrand;
	}

	/**
	 * The intron represented by this junction
	 * @return 
	 */
	shared_ptr<Intron> getIntron() const {
		return intron;
	}

	/**
	 * The intron size
	 * @return 
	 */
	int32_t getIntronSize() const {
		return intron != nullptr ? intron->size() : 0;
	}

	/**
	 * The start site of the left anchor represented by this junction
	 * @return 
	 */
	int32_t getLeftAncStart() const {
		return leftAncStart;
	}

	/**
	 * The end site of the right anchor represented by this junction
	 * @return 
	 */
	int32_t getRightAncEnd() const {
		return rightAncEnd;
	}

	/**
	 * The size of this junction from the start of the left exon anchor to the end
	 * of the right exon anchor
	 * @return 
	 */
	size_t size() const {
		return rightAncEnd - leftAncStart + 1;
	}

	/**
	 * The size of the left exon anchor
	 * @return 
	 */
	int32_t getLeftAnchorSize() const {
		return intron != nullptr ? intron->start - leftAncStart : 0;
	}

	/**
	 * The size of the right exon anchor
	 * @return 
	 */
	int32_t getRightAnchorSize() const {
		return intron != nullptr ? rightAncEnd - intron->end : 0;
	}

	/**
	 * The unique identifier assigned to this junction
	 * @return 
	 */
	uint32_t getId() const {
		return id;
	}

	/**
	 * Whether or not the user has marked this junction as genuine or not
	 * @return 
	 */
	bool isGenuine() const {
		return genuine;
	}

	/**
	 * The score assigned to this junction.  This could either be user-defined
	 * or set by the filtering tool.
	 * @return 
	 */
	double getScore() const {
		return score;
	}

	/**
	 * Whether or not there is a canonical donor and acceptor motif at the two base
	 * pairs at the start and end of the junction / intron
	 * @return 
	 */
	bool hasCanonicalSpliceSites() const {
		return this->canonicalSpliceSites == CanonicalSS::CANONICAL;
	}

	/**
	 * The type of junction this is: canonical, semi-canonical or non-canonical
	 * @return 
	 */
	CanonicalSS getSpliceSiteType() const {
		return this->canonicalSpliceSites;
	}

	/**
	 * The total number of spliced alignments directly supporting this junction
	 * @return 
	 */
	uint32_t getNbSplicedAlignments() const {
		return this->nbAlRaw;
	}

	/**
	 * The number of distinct alignments supporting this junction
	 * @return 
	 */
	uint32_t getNbDistinctAlignments() const {
		return this->nbAlDistinct;
	}

	/**
	 * The number of alignments convering only a single intron
	 * @return 
	 */
	uint32_t getNbUniquelySplicedAlignments() const {
		return this->nbAlRaw - this->nbAlMultiplySpliced;
	}

	/**
	 * The number of alignments covering 2 or more introns
	 * @return 
	 */
	uint32_t getNbMultiplySplicedAlignments() const {
		return this->nbAlMultiplySpliced;
	}

	/**
	 * The number of reads that probably align uniquely to the genome
	 * @return 
	 */
	uint32_t getNbUniquelyMappedAlignments() const {
		return this->nbAlUniquelyMapped;
	}

	/**
	 * The number of reads that probably align to multiple sites on the genome
	 * @return 
	 */
	uint32_t getNbMultiplyMappedAlignments() const {
		return this->nbAlRaw - this->nbAlUniquelyMapped;
	}

	/**
	 * The number of alignments that have a properly aligned pair as signalled
	 * by the properly paired SAM flag.  If running on single end data this will 
	 * always be 0.
	 * @return 
	 */
	uint32_t getNbBamProperlyPairedAlignments() const {
		return this->nbAlBamProperlyPaired;
	}

	/**
	 * The number of alignments that have a properly aligned pair as determined
	 * by portcullis.  If running on single end data this will always be 0.
	 * @return 
	 */
	uint32_t getNbPortcullisProperlyPairedAlignments() const {
		return this->nbAlPortcullisProperlyPaired;
	}

	/**
	 * The number of reliable split alignments supporting this junction.  With
	 * paired end data this will include split alignments that are both uniquely 
	 * mapping as well as properly paired.  For single end data this will be the 
	 * same as the number of uniquely mapping split alignments
	 * @return 
	 */
	uint32_t getNbReliableAlignments() const {
		return this->nbAlReliable;
	}

	/**
	 * The ratio of reliable alignments to raw alignments that support this junction
	 * @return 
	 */
	double getReliable2RawAlignmentRatio() const {
		return (double) this->nbAlReliable / (double) this->nbAlRaw;
	}

	/**
	 * The Shannon Entropy of this junction.  This is a measure of how well distributed
	 * the reads are around the junction.  The closely the alignment distribution
	 * is to a uniform distribution the higher the entropy score.
	 * @return 
	 */
	double getEntropy() const {
		return this->entropy;
	}

	/**
	 * The average number of mismatches of supporting alignments in this junction.
	 * @return 
	 */
	double getMeanMismatches() const {
		return this->meanMismatches;
	}

	/**
	 * The average read length of supporting alignments in this junction.
	 * @return 
	 */
	uint32_t getMeanReadLength() const {
		return this->meanReadLength;
	}

	/**
	 * The maximum of the shortest anchors from each supporting alignment
	 * @return
	 */
	uint32_t getMaxMinAnchor() const {
		return this->maxMinAnchor;
	}

	/**
	 * The maximum of the minimum matches on either side of the exon junction.
	 * This is similar to maxMinAnchor, except that it also includes mismatches
	 * into the calculation.
	 * @return 
	 */
	uint32_t getMaxMMES() const {
		return this->maxMMES;
	}

	/**
	 * The intron score is only calculated after filtering and represents the liklihood that
	 * the junction is invalid.  A score of 0 means the intron is around expected length.
	 * Values greater than 0 indicate the intron exceeds expected length.
	 * @return 
	 */
	double getIntronScore() const {
		return intronScore;
	}

	/**
	 * Hamming distance at the 3' end
	 * @return 
	 */
	int16_t getHammingDistance3p() const {
		return this->hammingDistance3p;
	}

	/**
	 * Hamming distance at the 5' end
	 * @return 
	 */
	int16_t getHammingDistance5p() const {
		return this->hammingDistance5p;
	}

	/**
	 * If true this junction shares no splice sites with other junctions in this
	 * dataset
	 * @return 
	 */
	bool isUniqueJunction() const {
		return this->uniqueJunction;
	}

	/**
	 * If true this junction is either unique (i.e. shares no splice sites with
	 * other junctions) or has the highest expression compared to other junctions
	 * that it shares splice sites with
	 * @return 
	 */
	bool isPrimaryJunction() const {
		return this->primaryJunction;
	}

	/**
	 * The number of junctions that can be directly traced downstream through spliced
	 * alignments
	 * @return 
	 */
	uint16_t getNbDownstreamJunctions() const {
		return nbDownstreamJunctions;
	}

	/**
	 * The number of junctions that can be directly traced upstream through spliced
	 * alignments
	 * @return 
	 */
	uint16_t getNbUpstreamJunctions() const {
		return nbUpstreamJunctions;
	}

	/**
	 * The distance to the next downstream junction that can be directly traced via 
	 * split alignments, or 0 if there no junctions directly downstream or if the
	 * next downstream junction is located within this junction
	 * @return 
	 */
	uint32_t getDistanceToNextDownstreamJunction() const {
		return distanceToNextDownstreamJunction;
	}

	/**
	 * The distance to the next upstream junction that can be directly traced via 
	 * split alignments, or 0 if there no junctions directly upstream or if the
	 * next upstream junction is located within this junction
	 * @return 
	 */
	uint32_t getDistanceToNextUpstreamJunction() const {
		return distanceToNextUpstreamJunction;
	}

	/**
	 * The distance to the nearest junction
	 * @return 
	 */
	uint32_t getDistanceToNearestJunction() const {
		return distanceToNearestJunction;
	}

	/**
	 * The Multiple mapping score (reflects mapping ambiguity.  Small score
	 * implies reads in this junction could align to other splice junctions across the
	 * genome).  This will be 0 unless the user requested extra processing.
	 * @return 
	 */
	double getMultipleMappingScore() const {
		return multipleMappingScore;
	}

	/**
	 * A score reflecting the coverage of flanking alignments around the junction.  Uses unspliced reads
	 * around both donor and acceptor sites.  If this is geniune junction you expect
	 * unspliced reads to drop off closer to the donor and acceptor sites you get
	 * (relative to read length).  We compare a window one read away from the d/a
	 * site and one read up to the d/a site.  We then substract the mean coverage from
	 * the more distant window to the close window.  We then add these score from
	 * both donor and acceptor sites together.  For genuine junctions this should be
	 * the score should be relatively large.  See TrueSight paper for more info.
	 * Only used if the user requested extra processing.
	 * @return 
	 */
	double getCoverage() const {
		return coverage;
	}

	/**
	 * The number of upstream non-spliced supporting reads.  Only used if the 
	 * user requested extra processing.
	 * @return 
	 */
	uint32_t getNbUpstreamFlankingAlignments() const {
		return nbUpstreamFlankingAlignments;
	}

	/**
	 * The number of downstream non-spliced supporting reads.  Only used if the
	 * user requested extra processing.
	 * @return 
	 */
	uint32_t getNbDownstreamFlankingAlignments() const {
		return nbDownstreamFlankingAlignments;
	}

	/**
	 * Whether or not this junction looks suspicious (i.e. it may not be genuine)
	 * due to no anchors extending beyond the first mismatch and if that location
	 * in either anchor is within 20bp of the junction.
	 * @return 
	 */
	bool isSuspicious() const {
		return this->suspicious;
	}

	/**
	 * Whether or not this junction is a potential false positive (i.e. may not be
	 * genuine).  For this to be true the junction must first be marked as suspicious
	 * then it must also be more than 99% likely that the maxmmes should have exceeded
	 * a given length, given the number of junctions supporting the junction but didn't.
	 * @return 
	 */
	bool isPotentialFalsePositive() const {
		return this->pfp;
	}

	/**
	 * Gets the junction anchor depth at the given distance from the intron
	 * @param index
	 * @return 
	 */
	uint32_t getJunctionAnchorDepth(size_t index) const {
		return junctionAnchorDepth[index];
	}


	/**
	 * This calls the relevant getter for the given name
	 * @param name
	 * @return 
	 */
	double getValueFromName(const string& name) const;



	// ***** Setters *****
	// Not all variables have setters

	void setId(uint32_t id) {
		this->id = id;
	}

	void setGenuine(bool genuine) {
		this->genuine = genuine;
	}

	void setScore(double score) {
		this->score = score;
	}

	void setDa1(string da1) {
		this->da1 = da1;
	}

	void setDa2(string da2) {
		this->da2 = da2;
	}
	
	void setNbSplicedAlignments(uint32_t nbAlRaw) {
		this->nbAlRaw = nbAlRaw;
	}

	void setNbDistinctAlignments(uint32_t nbDistinctAlignments) {
		this->nbAlDistinct = nbDistinctAlignments;
	}

	void setNbMultiplySplicedAlignments(uint32_t nbMultipleSplicedReads) {
		this->nbAlMultiplySpliced = nbMultipleSplicedReads;
	}

	void setNbBamProperlyPairedAlignments(uint32_t nbProperlyPairedReads) {
		this->nbAlBamProperlyPaired = nbProperlyPairedReads;
	}

	void setNbPortcullisProperlyPairedAlignments(uint32_t nbProperlyPairedReads) {
		this->nbAlPortcullisProperlyPaired = nbProperlyPairedReads;
	}

	void setNbUniquelyMappedAlignments(uint32_t nbUniquelyMappedReads) {
		this->nbAlUniquelyMapped = nbUniquelyMappedReads;
	}

	void setNbReliableAlignments(uint32_t nbReliableAlignments) {
		this->nbAlReliable = nbReliableAlignments;
	}

	void setEntropy(double entropy) {
		this->entropy = entropy;
	}

	void setMeanMismatches(double meanMismatches) {
		this->meanMismatches = meanMismatches;
	}

	void setMeanReadLength(uint32_t meanQueryLength) {
		this->meanReadLength = meanQueryLength;
	}

	void setMaxMinAnchor(int32_t maxMinAnchor) {
		this->maxMinAnchor = maxMinAnchor;
	}

	void setMaxMMES(uint32_t maxMMES) {
		this->maxMMES = maxMMES;
	}

	void setIntronScore(double intronScore) {
		this->intronScore = intronScore;
	}

	void setHammingDistance3p(int16_t hammingDistance3p) {
		this->hammingDistance3p = hammingDistance3p;
	}

	void setHammingDistance5p(int16_t hammingDistance5p) {
		this->hammingDistance5p = hammingDistance5p;
	}

	void setUniqueJunction(bool uniqueJunction) {
		this->uniqueJunction = uniqueJunction;
	}

	void setPrimaryJunction(bool primaryJunction) {
		this->primaryJunction = primaryJunction;
	}

	void setNbDownstreamJunctions(uint16_t nbDownstreamJunctions) {
		this->nbDownstreamJunctions = nbDownstreamJunctions;
	}

	void setNbUpstreamJunctions(uint16_t nbUpstreamJunctions) {
		this->nbUpstreamJunctions = nbUpstreamJunctions;
	}

	void setDistanceToNearestJunction(uint32_t distanceToNearestJunction) {
		this->distanceToNearestJunction = distanceToNearestJunction;
	}

	void setDistanceToNextDownstreamJunction(uint32_t distanceToNextDownstreamJunction) {
		this->distanceToNextDownstreamJunction = distanceToNextDownstreamJunction;
	}

	void setDistanceToNextUpstreamJunction(uint32_t distanceToNextUpstreamJunction) {
		this->distanceToNextUpstreamJunction = distanceToNextUpstreamJunction;
	}

	void setMultipleMappingScore(double multipleMappingScore) {
		this->multipleMappingScore = multipleMappingScore;
	}

	void setCoverage(double coverage) {
		this->coverage = coverage;
	}

	void setNbDownstreamFlankingAlignments(uint32_t nbDownstreamFlankingAlignments) {
		this->nbDownstreamFlankingAlignments = nbDownstreamFlankingAlignments;
	}

	void setNbUpstreamFlankingAlignments(uint32_t nbUpstreamFlankingAlignments) {
		this->nbUpstreamFlankingAlignments = nbUpstreamFlankingAlignments;
	}

	void setSuspicious(bool suspicious) {
		this->suspicious = suspicious;
	}

	void setPotentialFalsePositive(bool pfp) {
		this->pfp = pfp;
	}

	void setTrimmedLogDevCov(vector<double> trimmedLogDevCov) {
		this->trimmedLogDevCov = trimmedLogDevCov;
	}

	void setJunctionAnchorDepth(size_t index, uint32_t val) {
		junctionAnchorDepth[index] = val;
	}



	// ****** Methods for building junction anchors ******

	/**
	 * Add an alignment to this junction and update any associated properties
	 * @param al
	 */
	void addJunctionAlignment(const BamAlignment& al);

	/**
	 * Sets the canonical status of this junction based on the dinucleotides at the donor and acceptor
	 * sites.  This should also update the strand properties accordingly.
	 * @param seq1
	 * @param seq2
	 * @return 
	 */
	CanonicalSS setDonorAndAcceptorMotif(string seq1, string seq2);

	/**
	 * Whether or not the provided junction shares any splice sites with this
	 * junction
	 * @param other
	 * @return 
	 */
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


	// ****** Methods intended to be executed after all alignments are added to the junction

	/**
	 * This will iterate over the alignments in this junction and determine if there's consistency between
	 * the strand from the alignments
	 */
	void determineStrandFromReads();

	/**
	 * Extracts genomic content around this junction and updates any associated properties
	 * @param genomeMapper
	 */
	void processJunctionWindow(const GenomeMapper& genomeMapper);

	/**
	 * Process non-spliced alignments within a certain region upstream and downstream of the junction
	 * @param reader
	 * @param refLength
	 * @param meanQueryLength
	 * @param maxQueryLength
	 */
	void processJunctionVicinity(BamReader& reader, int32_t refLength, int32_t meanQueryLength, int32_t maxQueryLength);

	/**
	 * Based on the alignments in this junction calculate junction metrics
	 */
	void calcMetrics();

	/**
	 * Based on the alignments in this junction calculate junction metrics.  Use the given orientation
	 * to calculate properly paired reads.
	 * @param orientation
	 */
	void calcMetrics(Orientation orientation);

	/**
	 * Shannon Entropy (definition from "Graveley et al, The developmental 
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
	void calcAlignmentStats(Orientation orientation);

	/**
	 * Metric 13 and 14: Calculates the 5' and 3' hamming distances from a genomic
	 * region represented by this junction
	 */
	void calcHammingScores(const string& leftAnchor, const string& leftIntron,
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


	double calcCoverage(int32_t a, int32_t b, const vector<uint32_t>& coverageLevels);

	double calcCoverage(const vector<uint32_t>& coverageLevels);

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


	/**
	 * Given the markov models for the exons and the introns, calculate the potential coding score for
	 * this junction
	 * @param gmap
	 * @param exon
	 * @param intron
	 * @return 
	 */
	double calcCodingPotential(GenomeMapper& gmap, KmerMarkovModel& exon, KmerMarkovModel& intron) const;

	SplicingScores calcSplicingScores(GenomeMapper& gmap, KmerMarkovModel& donorT, KmerMarkovModel& donorF,
			KmerMarkovModel& acceptorT, KmerMarkovModel& acceptorF,
			PosMarkovModel& donorP, PosMarkovModel& acceptorP) const;


	/**
	 * Calculate the log deviation for the junction anchor depth count at a given location
	 * @param i
	 * @return 
	 */
	double calcJunctionAnchorDepthLogDeviation(size_t i) const;





	// **** Output methods ****

	/**
	 * Represent this junctions location as a string
	 * @return 
	 */
	string locationAsString() const {
		return this->intron->toString() + strandToChar(this->consensusStrand);
	}

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
	void outputIntronGFF(std::ostream &strm, const string& source);

	/**
	 * Complete human readable description of this junction
	 * @param strm
	 */
	void outputJunctionGFF(std::ostream &strm, const string& source);

	/**
	 * Complete human readable description of this junction
	 * @param strm
	 */
	void outputBED(std::ostream &strm, const string& prefix, bool bedscore);

	/**
	 * Represents this junction as a table row
	 * @param seq1
	 * @param seq2
	 * @return 
	 */
	friend ostream& operator<<(ostream &strm, Junction& j) {
		strm << j.id << "\t"
				<< *(j.intron) << "\t"
				<< j.getIntronSize() << "\t"
				<< j.leftAncStart << "\t"
				<< j.rightAncEnd << "\t"

				<< strandToChar(j.readStrand) << "\t"
				<< strandToChar(j.ssStrand) << "\t"
				<< strandToChar(j.consensusStrand) << "\t"

				<< j.da1 << "\t"
				<< j.da2 << "\t"
				<< cssToChar(j.canonicalSpliceSites) << "\t"
				
				<< j.score << "\t"
				<< j.suspicious << "\t"
				<< j.pfp << "\t"

				<< j.nbAlRaw << "\t"
				<< j.nbAlDistinct << "\t"
				<< j.getNbUniquelySplicedAlignments() << "\t"
				<< j.nbAlMultiplySpliced << "\t"
				<< j.nbAlUniquelyMapped << "\t"
				<< j.getNbMultiplyMappedAlignments() << "\t"
				<< j.nbAlBamProperlyPaired << "\t"
				<< j.nbAlPortcullisProperlyPaired << "\t"
				<< j.nbAlReliable << "\t"
				<< j.getReliable2RawAlignmentRatio() << "\t"

				<< j.entropy << "\t"
				<< j.meanMismatches << "\t"
				<< j.meanReadLength << "\t"
				<< j.maxMinAnchor << "\t"
				<< j.maxMMES << "\t"
				<< j.intronScore << "\t"

				<< j.hammingDistance5p << "\t"
				<< j.hammingDistance3p << "\t"

				<< j.uniqueJunction << "\t"
				<< j.primaryJunction << "\t"
				<< j.nbUpstreamJunctions << "\t"
				<< j.nbDownstreamJunctions << "\t"
				<< j.distanceToNextUpstreamJunction << "\t"
				<< j.distanceToNextDownstreamJunction << "\t"
				<< j.distanceToNearestJunction << "\t"

				<< j.multipleMappingScore << "\t"
				<< j.coverage << "\t"
				<< j.nbUpstreamFlankingAlignments << "\t"
				<< j.nbDownstreamFlankingAlignments;
				
		for (size_t i = 0; i < JAD_NAMES.size(); i++) {
			strm << "\t" << j.junctionAnchorDepth[i];
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
	
	static const vector<string> METRIC_NAMES;
	static const vector<string> JAD_NAMES;
	static const vector<string> STRAND_NAMES;	
};


const vector<string> Junction::METRIC_NAMES({
	"canonical_ss",
	"nb_raw_aln",
	"nb_dist_aln",
	"nb_us_aln",
	"nb_ms_aln",
	"nb_um_aln",
	"nb_mm_aln",
	"nb_bpp_aln",
	"nb_ppp_aln",
	"nb_rel_aln",
	"rel2raw",
	"entropy",
	"mean_mismatches",
	"mean_readlen",
	"max_min_anc",
	"maxmmes",
	"intron_score",
	"hamming5p",
	"hamming3p",
	"uniq_junc",
	"primary_junc",
	"nb_up_juncs",
	"nb_down_juncs",
	"dist_2_up_junc",
	"dist_2_down_junc",
	"dist_nearest_junc"
	"mm_score",
	"coverage",
	"up_aln",
	"down_aln"
});

typedef double (Junction::*junction_method_t)() const;
typedef std::map<std::string, junction_method_t> JuncFuncMap;

const JuncFuncMap JunctionFunctionMap = {
	{"nb_raw_aln", &Junction::getNbSplicedAlignments},
	{"nb_dist_aln", &Junction::getNbDistinctAlignments}
	/*"nb_us_aln",
	"nb_ms_aln",
	"nb_um_aln",
	"nb_mm_aln",
	"nb_bpp_aln",
	"nb_ppp_aln",
	"nb_rel_aln",
	"rel2raw",
	"entropy",
	"mean_mismatches",
	"mean_readlen",
	"max_min_anc",
	"maxmmes",
	"intron_score",
	"hamming5p",
	"hamming3p",
	"uniq_junc",
	"primary_junc",
	"nb_up_juncs",
	"nb_down_juncs",
	"dist_2_up_junc",
	"dist_2_down_junc",
	"dist_nearest_junc"
	"mm_score",
	"coverage",
	"up_aln",
	"down_aln"	*/
};

const vector<string> Junction::JAD_NAMES({
	"JAD01",
	"JAD02",
	"JAD03",
	"JAD04",
	"JAD05",
	"JAD06",
	"JAD07",
	"JAD08",
	"JAD09",
	"JAD10",
	"JAD11",
	"JAD12",
	"JAD13",
	"JAD14",
	"JAD15",
	"JAD16",
	"JAD17",
	"JAD18",
	"JAD19",
	"JAD20"
});

const vector<string> STRAND_NAMES = {
	"read-strand",
	"ss-strand",
	"consensus-strand"
};

struct JunctionComparator {

	inline bool operator()(const shared_ptr<Junction>& j1, const shared_ptr<Junction>& j2) const {
		return IntronComparator()(*(j1->getIntron()), *(j2->getIntron()));
	}
};

typedef shared_ptr<Junction> JunctionPtr;
typedef std::vector<shared_ptr<Junction> > JunctionList;

}
