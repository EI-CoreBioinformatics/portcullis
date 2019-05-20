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

#include <portcullis/junction.hpp>

const vector<string> portcullis::Junction::METRIC_NAMES({
	"canonical_ss",
	"score",
	"suspicious",
	"pfp",
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
	"nb_r1_pos",
	"nb_r1_neg",
	"nb_r2_pos",
	"nb_r2_neg",
	"entropy",
	"mean_mismatches",
	"mean_readlen",
	"max_min_anc",
	"maxmmes",
	"intron_score",
	"hamming5p",
	"hamming3p",
	"coding",
	"pws",
	"splice_sig",
	"uniq_junc",
	"primary_junc",
	"nb_up_juncs",
	"nb_down_juncs",
	"dist_2_up_junc",
	"dist_2_down_junc",
	"dist_nearest_junc",
	"mm_score",
	"coverage",
	"up_aln",
	"down_aln",
    "nb_samples"
});

const vector<string> portcullis::Junction::JAD_NAMES({
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

const vector<string> portcullis::Junction::STRAND_NAMES = {
	"read-strand",
	"ss-strand",
	"consensus-strand"
};

void portcullis::AlignmentInfo::calcMatchStats(const Intron& i, const uint32_t leftStart, const uint32_t rightEnd, const string& ancLeft, const string& ancRight) {

    if (leftStart > std::numeric_limits<int32_t>::max()) {
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Left start value is too large to process for junction: ") + i.toString()));
    }

    if (rightEnd > std::numeric_limits<int32_t>::max()) {
        BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Right end value is too large to process for junction: ") + i.toString()));
    }


    int32_t leftEnd = i.start - 1;
	int32_t rightStart = i.end + 1;
	int32_t qLeftStart = (int32_t)leftStart;
	int32_t qLeftEnd = leftEnd;
	int32_t qRightStart = rightStart;
	int32_t qRightEnd = (int32_t)rightEnd;

    string query = ba->getQuerySeq();
    if (query.size() <= 1) {
        // In this case the genome and query sequences do not correspond with one another.  Most
        // likely the cause of this is that the query sequence is not present in the alignment.
        // In which case just assume everything is fine.
		totalUpstreamMismatches = 0;
        totalDownstreamMismatches = 0;
        totalUpstreamMatches = leftEnd - leftStart + 1;
        totalDownstreamMatches = rightEnd - rightStart + 1;
        nbMismatches = totalUpstreamMismatches + totalDownstreamMismatches;
        upstreamMatches = totalUpstreamMismatches;
        downstreamMatches = totalDownstreamMismatches;
        minMatch = min(upstreamMatches, downstreamMatches);
        maxMatch = max(upstreamMatches, downstreamMatches);
        mmes = min(totalUpstreamMatches, totalDownstreamMatches);
    }
    else {

        string qAnchorLeft = ba->getPaddedQuerySeq(query, leftStart, leftEnd, qLeftStart, qLeftEnd, false);
        string qAnchorRight = ba->getPaddedQuerySeq(query, rightStart, rightEnd, qRightStart, qRightEnd, false);
        string gAnchorLeft = ba->getPaddedGenomeSeq(ancLeft, leftStart, leftEnd, qLeftStart, qLeftEnd, false);
        string gAnchorRight = ba->getPaddedGenomeSeq(ancRight, rightStart, rightEnd, qRightStart, qRightEnd, false);
        bool error = false;
        if (qAnchorLeft.size() != gAnchorLeft.size() || qAnchorLeft.empty()) {
            error = true;
            std::cerr << endl << "WARNING:  Skipping problematic alignment:" << endl
                                  <<  "Left anchor region for query and genome are not the same size." << endl
                                  << "Intron: " + i.toString() << endl
                                  << "Junction anchor limits: " + lexical_cast<string>(leftStart) + "," + lexical_cast<string>(rightEnd) << endl
                                  << "Genomic sequence: " + ancLeft << endl
                                  << "Alignment coords (before soft clipping): " + ba->toString(false) << endl
                                  << "Alignment coords (after soft clipping): " + ba->toString(true) << endl
                                  << "Read name: " + ba->deriveName() << endl
                                  << "Read seq (before soft clipping): " + ba->getQuerySeq() + " (" + lexical_cast<string>(ba->getQuerySeq().size()) + ")" << endl
                                  << "Read seq (after soft clipping): " + ba->getQuerySeqAfterClipping() + " (" + lexical_cast<string>(ba->getQuerySeqAfterClipping().size()) + ")" << endl
                                  << "Cigar: " + ba->getCigarAsString() << endl
                                  << "Left Anchor query seq:  \n" + qAnchorLeft + " (" + lexical_cast<string>(qAnchorLeft.size()) + ")" << endl
                                  << "Left Anchor genome seq: \n" + gAnchorLeft + " (" + lexical_cast<string>(gAnchorLeft.size()) + ")" << endl << endl;
        }
        if (qAnchorRight.size() != gAnchorRight.size() || qAnchorRight.empty()) {
            error = true;
            std::cerr << endl << "WARNING:  Skipping problematic alignment:" << endl
                                  << "Right Anchor region for query and genome are not the same size." << endl
                                  << "Intron: " + i.toString() << endl
                                  << "Junction anchor limits: " + lexical_cast<string>(leftStart) + "," + lexical_cast<string>(rightEnd) << endl
                                  << "Genomic sequence: " + ancRight << endl
                                  << "Alignment coords (before soft clipping): " + ba->toString(false) << endl
                                  << "Alignment coords (after soft clipping): " + ba->toString(true) << endl
                                  << "Read name: " + ba->deriveName() << endl
                                  << "Read seq (before soft clipping): " + ba->getQuerySeq() + " (" + lexical_cast<string>(ba->getQuerySeq().size()) + ")" << endl
                                  << "Read seq (after soft clipping): " + ba->getQuerySeqAfterClipping() + " (" + lexical_cast<string>(ba->getQuerySeqAfterClipping().size()) + ")" << endl
                                  << "Cigar: " + ba->getCigarAsString() << endl
                                  << "Right Anchor query seq:  \n" + qAnchorRight + " (" + lexical_cast<string>(qAnchorRight.size()) + ")" << endl
                                  << "Right Anchor genome seq: \n" + gAnchorRight + " (" + lexical_cast<string>(gAnchorRight.size()) + ")" << endl << endl;
        }
        if (!error) {
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
    }
}

uint32_t portcullis::AlignmentInfo::getNbMatchesFromStart(const string& query, const string& anchor) {
	for (size_t i = 0; i < query.size(); i++) {
		if (query[i] != anchor[i]) {
			return i;
		}
	}
	return query.size();
}

uint32_t portcullis::AlignmentInfo::getNbMatchesFromEnd(const string& query, const string& anchor) {
	for (size_t j = query.size(); j > 0; j--) {
		size_t i = j-1;
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
		return CanonicalSS::CANONICAL;
	}
	else if (seq == SEMI_CANONICAL_SEQ1 || seq == SEMI_CANONICAL_SEQ1_RC ||
			 seq == SEMI_CANONICAL_SEQ2 || seq == SEMI_CANONICAL_SEQ2_RC) {
		return CanonicalSS::SEMI_CANONICAL;
	}
	else {
		return CanonicalSS::NO;
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
	}// Check for these after canonicals for performance reasons
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
	id = 0;
	leftAncStart = _leftAncStart;
	rightAncEnd = _rightAncEnd;
	readStrand = Strand::UNKNOWN;
	ssStrand = Strand::UNKNOWN;
	consensusStrand = Strand::UNKNOWN;
	genuine = false;
	score = 0.0;
	suspicious = false;
	pfp = false;
	canonicalSpliceSites = CanonicalSS::NO;
	nbAlRaw = 0;
	nbAlDistinct = 0;
	nbAlMultiplySpliced = 0;
	nbAlUniquelyMapped = 0;
	nbAlBamProperlyPaired = 0;
	nbAlPortcullisProperlyPaired = 0;
	nbAlReliable = 0;
	nbAlR1Pos = 0;
	nbAlR1Neg = 0;
	nbAlR2Pos = 0;
	nbAlR2Neg = 0;
	entropy = 0;
	meanMismatches = 0;
	meanReadLength = 0;
	maxMinAnchor = intron->minAnchorLength(_leftAncStart, _rightAncEnd);
	maxMMES = 0;
	intronScore = 0.0;
	hammingDistance5p = 10;
	hammingDistance3p = 10;
	codingPotential = 0.0;
	positionWeightScore = 0.0;
	splicingSignal = 0.0;
	uniqueJunction = false;
	primaryJunction = false;
	nbUpstreamJunctions = 0;
	nbDownstreamJunctions = 0;
	distanceToNextUpstreamJunction = 0;
	distanceToNextDownstreamJunction = 0;
	distanceToNearestJunction = 0;
	coverage = 0.0;
	multipleMappingScore = 0.0;
	nbUpstreamFlankingAlignments = 0;
	nbDownstreamFlankingAlignments = 0;
    nbSamples = 1;
	junctionAnchorDepth.clear();
	for (size_t i = 0; i < JAD_NAMES.size(); i++) {
		junctionAnchorDepth.push_back(0);
	}
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
	readStrand = j.readStrand;
	ssStrand = j.ssStrand;
	consensusStrand = j.consensusStrand;
	id = j.id;
	suspicious = j.suspicious;
	pfp = j.pfp;
	genuine = j.genuine;
	score = j.score;
	canonicalSpliceSites = j.canonicalSpliceSites;
	nbAlRaw = j.nbAlRaw;
	nbAlDistinct = j.nbAlDistinct;
	nbAlMultiplySpliced = j.nbAlMultiplySpliced;
	nbAlUniquelyMapped = j.nbAlUniquelyMapped;
	nbAlBamProperlyPaired = j.nbAlBamProperlyPaired;
	nbAlPortcullisProperlyPaired = j.nbAlPortcullisProperlyPaired;
	nbAlReliable = j.nbAlReliable;
	nbAlR1Pos = j.nbAlR1Pos;
	nbAlR1Neg = j.nbAlR1Neg;
	nbAlR2Pos = j.nbAlR2Pos;
	nbAlR2Neg = j.nbAlR2Neg;
	entropy = j.entropy;
	meanMismatches = j.meanMismatches;
	meanReadLength = j.meanReadLength;
	maxMinAnchor = j.maxMinAnchor;
	maxMMES = j.maxMMES;
	intronScore = j.intronScore;
	hammingDistance5p = j.hammingDistance5p;
	hammingDistance3p = j.hammingDistance3p;
	codingPotential = j.codingPotential;
	positionWeightScore = j.positionWeightScore;
	splicingSignal = j.splicingSignal;
	uniqueJunction = j.uniqueJunction;
	primaryJunction = j.primaryJunction;
	nbUpstreamJunctions = j.nbUpstreamJunctions;
	nbDownstreamJunctions = j.nbDownstreamJunctions;
	distanceToNextUpstreamJunction = j.distanceToNextUpstreamJunction;
	distanceToNextDownstreamJunction = j.distanceToNextDownstreamJunction;
	distanceToNearestJunction = j.distanceToNearestJunction;
	coverage = j.coverage;
	multipleMappingScore = j.multipleMappingScore;
	nbUpstreamFlankingAlignments = j.nbUpstreamFlankingAlignments;
	nbDownstreamFlankingAlignments = j.nbDownstreamFlankingAlignments;
    nbSamples = j.nbSamples;
	if (withAlignments) {
		for (size_t i = 0; i < j.alignments.size(); i++) {
			this->alignments.push_back(make_shared<AlignmentInfo>(j.alignments[i]->ba));
			this->alignmentCodes.push_back(j.alignments[i]->nameCode);
		}
	}
	trimmedCoverage.clear();
	for (auto & x : j.trimmedCoverage) {
		trimmedCoverage.push_back(x);
	}
	trimmedLogDevCov.clear();
	for (auto & x : j.trimmedLogDevCov) {
		trimmedLogDevCov.push_back(x);
	}
	junctionAnchorDepth.clear();
	for (size_t i = 0; i < JAD_NAMES.size(); i++) {
		junctionAnchorDepth.push_back(j.getJunctionAnchorDepth(i));
	}
}

// **** Destructor ****

portcullis::Junction::~Junction() {
	alignments.clear();
	junctionAnchorDepth.clear();
}

void portcullis::Junction::clearAlignments() {
	alignments.clear();
}

void portcullis::Junction::addJunctionAlignment(const BamAlignment& al) {
	// Make sure we take a proper copy of this alignment for safe storage
	AlignmentInfoPtr aip = make_shared<AlignmentInfo>(make_shared<BamAlignment>(al));
	this->alignments.push_back(aip);
	this->alignmentCodes.push_back(aip->nameCode);
	this->nbAlRaw = this->alignments.size();
	if (al.isFirstMate()) {
		if (!al.isReverseStrand()) {
			this->nbAlR1Pos++;
		}
		else {
			this->nbAlR1Neg++;
		}
	}
	else {
		if (!al.isReverseStrand()) {
			this->nbAlR2Pos++;
		}
		else {
			this->nbAlR2Neg++;
		}
	}
	if (al.getNbJunctionsInRead() > 1) {
		this->nbAlMultiplySpliced++;
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
	this->da1 = this->consensusStrand == Strand::NEGATIVE ? SeqUtils::reverseComplement(seq2) : seq1;
	this->da2 = this->consensusStrand == Strand::NEGATIVE ? SeqUtils::reverseComplement(seq1) : seq2;
	return this->canonicalSpliceSites;
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
	uint32_t otherMinAnchor = intron->minAnchorLength(otherStart, otherEnd);
	maxMinAnchor = max(maxMinAnchor, otherMinAnchor);
}

void portcullis::Junction::determineStrandFromReads() {
	uint32_t nb_pos = 0;
	uint32_t nb_neg = 0;
	uint32_t nb_unk = 0;
	for (const auto & a : alignments) {
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
	if ((double) nb_pos / (double) total >= threshold) {
		readStrand = Strand::POSITIVE;
	}
	else if ((double) nb_neg / (double) total >= threshold) {
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
	string donor = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->start, intron->start + 1);
	string acceptor = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->end - 1, intron->end);
	int donorLen = donor.length();
	int acceptorLen = acceptor.length();
	if (donorLen == 0)
		BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Can't find donor site (left side splice site) region for junction: ") + this->intron->toString()));
	if (donorLen != 2)
		BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Retrieved sequence for left side splice site of junction ") + this->intron->toString() + " is not the expected length" +
							  "\nRetrieved sequence Length: " + lexical_cast<string>(donorLen) +
							  "\nExpected sequence length: " + lexical_cast<string>(2) + "\n"));
	if (acceptorLen == 0)
		BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Can't find acceptor site (right side splice site) region for junction: ") + this->intron->toString()));
	if (acceptorLen != 2)
		BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Retrieved sequence for right side splice site of junction ") + this->intron->toString() + " is not the expected length" +
							  "\nRetrieved sequence Length: " + lexical_cast<string>(acceptorLen) +
							  "\nExpected sequence length: " + lexical_cast<string>(2) + "\n"));
	boost::to_upper(donor); // Removes any lowercase bases representing repeats
	boost::to_upper(acceptor); // Removes any lowercase bases representing repeats
	this->setDonorAndAcceptorMotif(donor, acceptor);
	// Just access the whole junction region
	string leftAnc = genomeMapper.fetchBases(intron->ref.name.c_str(), leftAncStart, intron->start - 1);
	string rightAnc = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->end + 1, rightAncEnd);
	string leftInt = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->start, intron->start + 9);
	string rightInt = genomeMapper.fetchBases(intron->ref.name.c_str(), intron->end - 9, intron->end);
	int leftAncLen = leftAnc.length();
	int leftIntLen = leftInt.length();
	int rightAncLen = rightAnc.length();
	int rightIntLen = rightInt.length();
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
	for (const auto & a : alignments) {
		a->calcMatchStats(*getIntron(), this->getLeftAncStart(), this->getRightAncEnd(), leftAnc, rightAnc);
	}
	// MaxMMES can now use info in alignments
	this->calcMismatchStats();
}

void portcullis::Junction::processJunctionVicinity(BamReader& reader, int32_t refLength, int32_t maxQueryLength) {
	int32_t refId = intron->ref.index;
	uint32_t nbLeftFlankingAlignments = 0, nbRightFlankingAlignments = 0;
	int32_t regionStart = leftAncStart - maxQueryLength - 1;
	regionStart = regionStart < 0 ? 0 : regionStart;
	int32_t regionEnd = rightAncEnd + maxQueryLength + 1;
	regionEnd = regionEnd >= refLength ? refLength - 1 : regionEnd;
	// Focus only on the (expanded... to be safe...) region of interest
	reader.setRegion(refId, regionStart, regionEnd);
	while (reader.next()) {
		const BamAlignment& ba = reader.current();
		int32_t pos = ba.getStart();
		//TODO: Should we consider strand specific reads differently here?
		// Look for left flanking alignments
		if (intron->start > pos &&
				leftAncStart <= ba.getEnd()) {
			nbLeftFlankingAlignments++;
		}
		// Look for right flanking alignments
		if (rightAncEnd >= pos &&
				intron->end < pos) {
			nbRightFlankingAlignments++;
		}
	}
	this->nbUpstreamFlankingAlignments = nbLeftFlankingAlignments;
	this->nbDownstreamFlankingAlignments = nbRightFlankingAlignments;
}

void portcullis::Junction::calcMetrics() {
	this->calcMetrics(Orientation::UNKNOWN);
}

void portcullis::Junction::calcMetrics(Orientation orientation) {
	determineStrandFromReads();
	calcEntropy();
	calcAlignmentStats(orientation);
}

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
 * @param junctionPositions start index positions of each junction alignment
 *
 * @return The entropy of this junction
 */
double portcullis::Junction::calcEntropy() {
	vector<int32_t> junctionPositions;
	for (const auto & a : alignments) {
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
	for (size_t i = 0; i < nbJunctionAlignments; i++) {
		int32_t pos = junctionPositions[i];
		readsAtOffset++;
		if (pos != lastOffset || i == nbJunctionAlignments - 1) {
			double pI = (double) readsAtOffset / (double) nbJunctionAlignments;
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
void portcullis::Junction::calcAlignmentStats(Orientation orientation) {
	int32_t lastStart = -1, lastEnd = -1;
	nbAlDistinct = 0;
	nbAlReliable = 0;
	nbUpstreamJunctions = 0;
	nbDownstreamJunctions = 0;
	const bool properPairedCheck = doProperPairCheck(orientation);
	//cout << junctionAlignments.size() << endl;
	for (const auto & a : alignments) {
		BamAlignmentPtr ba = a->ba;
		const int32_t start = ba->getStart();
		const int32_t end = ba->getEnd();
		if (start != lastStart || end != lastEnd) {
			nbAlDistinct++;
			lastStart = start;
			lastEnd = end;
		}
		bool reliable = true;
		if (ba->getMapQuality() >= MAP_QUALITY_THRESHOLD) {
			nbAlUniquelyMapped++;
		}
		else {
			reliable = false;
		}
		// Get properly paired BAM flag regardless
		if (ba->isProperPair()) {
			nbAlBamProperlyPaired++;
		}
		if (properPairedCheck) {
			bool pp = ba->calcIfProperPair(orientation);
			if (pp) {
				nbAlPortcullisProperlyPaired++;
			}
			else {
				reliable = false;
			}
		}
		if (reliable) {
			nbAlReliable++;
		}
		uint32_t upjuncs = 0;
		uint32_t downjuncs = 0;
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
void portcullis::Junction::calcHammingScores(const string& leftAnchor, const string& leftIntron,
		const string& rightIntron, const string& rightAnchor) {
	const int32_t leftDelta = leftAnchor.size() - rightIntron.size();
	const int32_t leftOffset = leftDelta <= 0 ? 0 : leftDelta;
	const uint32_t leftLen = min(leftAnchor.size(), rightIntron.size());
	const uint32_t rightLen = min(leftIntron.size(), rightAnchor.size());
	// TODO, might want to modify this logic later, but worst case is that the 5' and
	// 3' results are swapped
	Strand s = consensusStrand != Strand::UNKNOWN ? consensusStrand : Strand::UNKNOWN;
	const string la = leftAnchor.size() > leftLen ? leftAnchor.substr(leftOffset, leftLen) : leftAnchor;
	const string li = leftIntron.size() > rightLen ? leftIntron.substr(0, rightLen) : leftIntron;
	const string ri = rightIntron.size() > leftLen ? rightIntron.substr(leftOffset, leftLen) : rightIntron;
	const string ra = rightAnchor.size() > rightLen ? rightAnchor.substr(0, rightLen) : rightAnchor;
	string anchor5p;
	string intron5p;
	string intron3p;
	string anchor3p;
	if (s == Strand::NEGATIVE) {
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
	for (const auto & a : alignments) {
		// Update maxMMES for this alignment
		maxMMES = max(maxMMES, a->mmes);
		// Update total number of mismatches in this junction
		nbMismatches += a->nbMismatches;
		// Keep a record of the first mismatch detected
		if (a->minMatch > 0) {
			firstMismatch = min(firstMismatch, a->minMatch);
		}
		// Update junction overhang vector
		for (uint16_t i = 0; i < JAD_NAMES.size() && i < a->minMatch; i++) {
			junctionAnchorDepth[i]++;
		}
	}
	// Set mean mismatches across junction
	meanMismatches = (double) nbMismatches / (double) alignments.size();
	// Assuming we have some mismatches determine if this junction has no overhangs
	// extending beyond first mismatch.  If so determine if that distance is small
	// enough to consider the junction as suspicious
	if (nbMismatches > 0 && firstMismatch < 20) {
		bool found = false;
		for (const auto & a : alignments) {
			if (a->minMatch > firstMismatch) {
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
	for (const auto & a : alignmentCodes) {
		M += map[a]; // Number of multiple splitting patterns
	}
	this->multipleMappingScore = (double) N / (double) M;
}

double portcullis::Junction::calcCoverage(int32_t a, int32_t b, const vector<uint32_t>& coverageLevels) {
	double multiplier = 1.0 / (b - a);
	uint32_t readCount = 0;
	for (int32_t i = a; i <= b; i++) {
		// Don't do anything stupid!
		if (i >= 0 && i < (int32_t)coverageLevels.size()) {
			readCount += coverageLevels[i];
		}
	}
	return multiplier * (double) readCount;
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
		calcCoverage(donorStart, donorMid - 1, coverageLevels) -
		calcCoverage(donorMid, donorEnd, coverageLevels);
	double acceptorCoverage =
		calcCoverage(acceptorMid, acceptorEnd, coverageLevels) -
		calcCoverage(acceptorStart, acceptorMid - 1, coverageLevels);
	coverage = donorCoverage + acceptorCoverage;
	return coverage;
}

double portcullis::Junction::calcIntronScore(const uint32_t threshold) {
	this->setIntronScore((uint32_t)this->intron->size() <= threshold ? 0.0 : log(this->intron->size() - threshold));
	return this->intronScore;
}

double portcullis::Junction::getValueFromName(const string& name) const {
	double val = 0.0;
	JuncUint32FuncMap::const_iterator uif = JunctionUint32FunctionMap.find(name);
	if (uif == JunctionUint32FunctionMap.end()) {
		JuncDoubleFuncMap::const_iterator df = JunctionDoubleFunctionMap.find(name);
		if (df == JunctionDoubleFunctionMap.end()) {
			JuncBoolFuncMap::const_iterator bf = JunctionBoolFunctionMap.find(name);
			if (bf == JunctionBoolFunctionMap.end()) {
				BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
										  "Unrecognised junction property: ") + name));
			}
			else {
				val = (double) (this->*(bf->second))();
			}
		}
		else {
			val = (this->*(df->second))();
		}
	}
	else {
		val = (double) (this->*(uif->second))();
	}
	return val;
}

bool portcullis::Junction::isNumericType(const string& name) {
	JuncUint32FuncMap::const_iterator uif = JunctionUint32FunctionMap.find(name);
	if (uif == JunctionUint32FunctionMap.end()) {
		JuncDoubleFuncMap::const_iterator df = JunctionDoubleFunctionMap.find(name);
		if (df == JunctionDoubleFunctionMap.end()) {
			JuncBoolFuncMap::const_iterator bf = JunctionBoolFunctionMap.find(name);
			if (bf == JunctionBoolFunctionMap.end()) {
				return false;
			}
		}
	}
	return true;
}

string portcullis::Junction::getStringFromName(const string& name) const {
	JuncStringFuncMap::const_iterator sf = JunctionStringFunctionMap.find(name);
	if (sf == JunctionStringFunctionMap.end()) {
		BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Unrecognised junction property: ") + name));
	}
	return (this->*(sf->second))();
}

bool portcullis::Junction::isStringType(const string& name) {
	JuncStringFuncMap::const_iterator sf = JunctionStringFunctionMap.find(name);
	if (sf == JunctionStringFunctionMap.end()) {
		return false;
	}
	return true;
}

// **** Output methods ****

/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::outputDescription(std::ostream &strm, string delimiter) {
	strm << "*** Intron ***" << delimiter;
	if (intron != nullptr) {
		intron->outputDescription(strm, delimiter);
		strm << delimiter << "Intron Size: " << getIntronSize();
	}
	else {
		strm << "No location set";
	}
	strm << delimiter
		 << "*** Anchors ***" << delimiter
		 << "Anchor limits: (" << leftAncStart << ", " << rightAncEnd << ")" << delimiter
		 << "Anchor sizes: (" << getLeftAnchorSize() << ", " << getRightAnchorSize() << ")" << delimiter
		 << "*** Strand ***" << delimiter
		 << "Reads Strand: " << strandToString(readStrand) << delimiter
		 << "Splice Site Strand: " << strandToString(ssStrand) << delimiter
		 << "Consensus Strand: " << strandToString(consensusStrand) << delimiter
		 << "*** Confidence ***" << delimiter
		 << "Canonical?: " << boolalpha << cssToString(canonicalSpliceSites) << "; Sequences: (" << da1 << " " << da2 << ")" << delimiter
		 << "Filter score: " << score << delimiter
		 << "Suspicious? (no anchors extending beyond first mismatch): " << boolalpha << suspicious << delimiter
		 << "Potential False Positive? (Suspicious and MaxMMES should have been greater given junction depth): " << boolalpha << pfp << delimiter
		 << "*** Alignment counts ***" << delimiter
		 << "# Total Spliced Alignments: " << nbAlRaw << delimiter
		 << "# Distinct Alignments: " << nbAlDistinct << delimiter
		 << "# Uniquely Spliced Alignments: " << getNbUniquelySplicedAlignments() << delimiter
		 << "# Multiply Spliced Alignments: " << nbAlMultiplySpliced << delimiter
		 << "# Uniquely Mapped Alignments: " << nbAlUniquelyMapped << delimiter
		 << "# Multiply Mapped Alignments: " << getNbMultiplyMappedAlignments() << delimiter
		 << "# Properly paired (bam flag): " << nbAlBamProperlyPaired << delimiter
		 << "# Properly paired (portcullis): " << nbAlPortcullisProperlyPaired << delimiter
		 << "# Reliable (MapQ >=" << MAP_QUALITY_THRESHOLD << " + portcullis properly paired) Alignments: " << nbAlReliable << delimiter
		 << "# R1 (+" << nbAlR1Pos << ",-" << nbAlR1Neg << "); # R2 (+" << nbAlR2Pos << ",-" << nbAlR2Neg << ")" << delimiter
		 << "*** RNA seq derived Junction stats ***" << delimiter
		 << "Entropy: " << entropy << delimiter
		 << "Mean mismatches: " << meanMismatches << delimiter
		 << "Mean read length: " << meanReadLength << delimiter
		 << "MaxMinAnchor: " << maxMinAnchor << delimiter
		 << "MaxMMES: " << maxMMES << delimiter
		 << "Intron score: " << intronScore << delimiter
		 << "*** Genome derived Junction stats ***" << delimiter
		 << "Hamming Distance 5': " << hammingDistance5p << delimiter
		 << "Hamming Distance 3': " << hammingDistance3p << delimiter
		 << "*** Junction group properties ***" << delimiter
		 << "Unique Junction: " << boolalpha << uniqueJunction << delimiter
		 << "Primary Junction: " << boolalpha << primaryJunction << delimiter
		 << "# Upstream Junctions: " << nbUpstreamJunctions << delimiter
		 << "# Downstream Junctions: " << nbDownstreamJunctions << delimiter
		 << "Distance to next upstream junction: " << distanceToNextUpstreamJunction << delimiter
		 << "Distance to next downstream junction: " << distanceToNextDownstreamJunction << delimiter
		 << "Distance to nearest junction: " << distanceToNearestJunction << delimiter
		 << "*** Extra metrics ***" << delimiter
		 << "Multiple mapping score: " << multipleMappingScore << delimiter
		 << "Coverage: " << coverage << delimiter
		 << "# Upstream Non-Spliced Alignments: " << nbUpstreamFlankingAlignments << delimiter
		 << "# Downstream Non-Spliced Alignments: " << nbDownstreamFlankingAlignments << delimiter
         << "# Samples: " << nbSamples;
}

/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::condensedOutputDescription(std::ostream &strm, string delimiter) {
	strm << "Strand: " << strandToString(consensusStrand) << delimiter
		 << "Canonical?=" << cssToString(canonicalSpliceSites) << delimiter
		 << "Score=" << score << delimiter
		 << "NbAlignments=" << getNbSplicedAlignments() << delimiter
		 << "NbDistinct=" << nbAlDistinct << delimiter
		 << "NbReliable=" << nbAlReliable << delimiter
		 << "Entropy=" << entropy << delimiter
		 << "MaxMMES=" << maxMMES << delimiter
		 << "HammingDistance5=" << hammingDistance5p << delimiter
		 << "HammingDistance3=" << hammingDistance3p << delimiter
		 << "UniqueJunction=" << boolalpha << uniqueJunction << delimiter
		 << "PrimaryJunction=" << boolalpha << primaryJunction << delimiter;
}

/**
 * Complete human readable description of this intron (for augustus hints)
 * @param strm
 */
void portcullis::Junction::outputIntronGFF(std::ostream &strm, const string& source) {
	// Use intron strand if known, otherwise use the predicted strand,
	// if predicted strand is also unknown then use "." to indicated unstranded
	const char strand = consensusStrand == Strand::UNKNOWN ? '?' : strandToChar(consensusStrand);
	string juncId = string("junc_") + lexical_cast<string>(id);
	// Modify coordinates to 1-based end inclusive
	// Output junction parent
	strm << intron->ref.name << "\t"
		 << source << "\t" // source
		 << "intron" << "\t" // type (may change later)
		 << intron->start + 1 << "\t" // start
		 << intron->end + 1 << "\t" // end
		 << nbAlRaw << "\t" // No score for the moment
		 << strand << "\t" // strand
		 << "." << "\t" // Just put "." for the phase
		 // Removing this as it causes issues with PASA downstream
		 /**<< "Note=cov:" << nbJunctionAlignments
		                            << "|rel:" << this->nbReliableAlignments
		                            << "|ent:" << std::setprecision(4) << this->entropy << std::setprecision(9)
		                            << "|maxmmes:" << this->maxMMES
		                            << "|ham:" << min(this->hammingDistance3p, this->hammingDistance5p) << ";"  // Number of times it was seen**/
		 << "mult=" << nbAlRaw << ";" // Coverage for augustus
		 << "grp=" << juncId << ";" // ID for augustus
		 << "src=E"; // Source for augustus
	strm << endl;
}

/**
 * Complete human readable description of this junction
 * @param strm
 */
void portcullis::Junction::outputJunctionGFF(std::ostream &strm, const string& source) {
	// Use intron strand if known, otherwise use the predicted strand,
	// if predicted strand is also unknown then use "." to indicated unstranded
	const char strand = consensusStrand == Strand::UNKNOWN ? '?' : strandToChar(consensusStrand);
	string juncId = string("junc_") + lexical_cast<string>(id);
	// Modify coordinates to 1-based end inclusive
	// Output junction parent
	strm << intron->ref.name << "\t"
		 << source << "\t" // source
		 << "match" << "\t" // type (may change later)
		 << leftAncStart + 1 << "\t" // start
		 << rightAncEnd + 1 << "\t" // end
		 << "0.0" << "\t" // No score for the moment
		 << strand << "\t" // strand
		 << "." << "\t" // Just put "." for the phase
		 << "ID=" << juncId << ";" // ID of the intron
		 << "Name=" << juncId << ";" // ID of the intron
		 << "Note=cov:" << nbAlRaw
		 << "|rel:" << this->nbAlReliable
		 << "|ent:" << std::setprecision(4) << this->entropy << std::setprecision(9)
		 << "|maxmmes:" << this->maxMMES
		 << "|ham:" << min(this->hammingDistance3p, this->hammingDistance5p) << ";" // Number of times it was seen
		 << "mult=" << nbAlRaw << ";" // Coverage for augustus
		 << "grp=" << juncId << ";" // ID for augustus
		 << "src=E;"; // Source for augustus
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
void portcullis::Junction::outputBED(std::ostream &strm, const string& prefix, bool bedscore) {
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
	strm << std::fixed << std::setprecision(3);
	// Output junction parent
	strm << intron->ref.name << "\t" // chrom
		 << leftAncStart << "\t" // chromstart
		 << rightAncEnd + 1 << "\t" // chromend (adding 1 as end position is exclusive)
		 << juncId << "\t" // name
		 << (bedscore ? this->getScore() : this->getNbSplicedAlignments()) << "\t" // Use the depth as the score for the moment
		 << strand << "\t" // strand
		 << intron->start << "\t" // thickstart
		 << intron->end + 1 << "\t" // thickend  (adding 1 as end position is exclusive)
		 << "255,0,0" << "\t" // Just use red for the moment
		 << "2" << "\t" // 2 blocks: Left and right block
		 << blockSizes << "\t"
		 << blockStarts << endl;
}



// **** Static methods ****

/**
 * Header for table output
 * @return
 */
string portcullis::Junction::junctionOutputHeader() {
	return string("index\t") + Intron::locationOutputHeader() + "\tsize\tleft\tright\t" +
		   boost::algorithm::join(Junction::STRAND_NAMES, "\t") + "\tss1\tss2\t" +
		   boost::algorithm::join(Junction::METRIC_NAMES, "\t") + "\t" +
		   boost::algorithm::join(Junction::JAD_NAMES, "\t");
}

shared_ptr<portcullis::Junction> portcullis::Junction::parse(const string& line) {
	vector<string> parts; // #2: Search for tokens
	boost::split(parts, line, boost::is_any_of("\t"), boost::token_compress_on);
	uint32_t expected_cols = 11 + Junction::STRAND_NAMES.size() + Junction::METRIC_NAMES.size() + Junction::JAD_NAMES.size();
	if (parts.size() != expected_cols) {
	  std::string sparts;
	  sparts = accumulate(begin(parts), end(parts), sparts);
		BOOST_THROW_EXCEPTION(JunctionException() << JunctionErrorInfo(string(
								  "Could not parse line due to incorrect number of columns.  This is probably a version mismatch.  Check file and portcullis versions.  Expected ")
							  + std::to_string(expected_cols) + " columns.  Found "
									       + std::to_string(parts.size()) + ".  Line:\n" + sparts));
	}
	// Create intron
	IntronPtr intron = make_shared<Intron>(
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
								 intron,
								 lexical_cast<int32_t>(parts[7]),
								 lexical_cast<int32_t>(parts[8])
							 );
	j->setId(lexical_cast<uint32_t>(parts[0]));
	// Index... saves having to remember the column numbers
	uint16_t i = 9;
	// Set predictions to junction
	j->readStrand = strandFromChar(parts[i++][0]);
	j->ssStrand = strandFromChar(parts[i++][0]);
	j->consensusStrand = strandFromChar(parts[i++][0]);
	// Splice site properties
	j->setDa1(parts[i++]);
	j->setDa2(parts[i++]);
	j->canonicalSpliceSites = cssFromChar(parts[i++][0]);
	// Confidence properties
	j->setScore(lexical_cast<double>(parts[i++]));
	j->setSuspicious(lexical_cast<bool>(parts[i++]));
	j->setPotentialFalsePositive(lexical_cast<bool>(parts[i++]));
	// Alignment counts
	j->setNbSplicedAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbDistinctAlignments(lexical_cast<uint32_t>(parts[i++]));
	i++; // uniquely spliced alignments not required
	j->setNbMultiplySplicedAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbUniquelyMappedAlignments(lexical_cast<uint32_t>(parts[i++]));
	i++; // multiply mapped alignments not required
	j->setNbBamProperlyPairedAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbPortcullisProperlyPairedAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbReliableAlignments(lexical_cast<uint32_t>(parts[i++]));
	i++; // reliable2raw ratio not required
	j->setNbR1PosAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbR1NegAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbR2PosAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbR2NegAlignments(lexical_cast<uint32_t>(parts[i++]));

	// RNAseq derived Junction stats
	j->setEntropy(lexical_cast<double>(parts[i++]));
	j->setMeanMismatches(lexical_cast<double>(parts[i++]));
	j->setMeanReadLength(lexical_cast<double>(parts[i++]));
	j->setMaxMinAnchor(lexical_cast<int32_t>(parts[i++]));
	j->setMaxMMES(lexical_cast<uint32_t>(parts[i++]));
	j->setIntronScore(lexical_cast<double>(parts[i++]));
	// Genome derived junction stats
	j->setHammingDistance5p(lexical_cast<uint32_t>(parts[i++]));
	j->setHammingDistance3p(lexical_cast<uint32_t>(parts[i++]));
	j->setCodingPotential(lexical_cast<double>(parts[i++]));
	j->setPositionWeightScore(lexical_cast<double>(parts[i++]));
	j->setSplicingSignal(lexical_cast<double>(parts[i++]));
	// Junction group properties
	j->setUniqueJunction(lexical_cast<bool>(parts[i++]));
	j->setPrimaryJunction(lexical_cast<bool>(parts[i++]));
	j->setNbUpstreamJunctions(lexical_cast<uint16_t>(parts[i++]));
	j->setNbDownstreamJunctions(lexical_cast<uint16_t>(parts[i++]));
	j->setDistanceToNextUpstreamJunction(lexical_cast<uint32_t>(parts[i++]));
	j->setDistanceToNextDownstreamJunction(lexical_cast<uint32_t>(parts[i++]));
	j->setDistanceToNearestJunction(lexical_cast<uint32_t>(parts[i++]));
	// Extra metrics requiring additional processing
	j->setMultipleMappingScore(lexical_cast<double>(parts[i++]));
	j->setCoverage(lexical_cast<double>(parts[i++]));
	j->setNbUpstreamFlankingAlignments(lexical_cast<uint32_t>(parts[i++]));
	j->setNbDownstreamFlankingAlignments(lexical_cast<uint32_t>(parts[i++]));
    j->setNbSamples(lexical_cast<uint32_t>(parts[i++]));
	// Read Junction anchor depths
	for (size_t k = 0; k < Junction::JAD_NAMES.size(); k++) {
		j->setJunctionAnchorDepth(k, lexical_cast<uint32_t>(parts[i + k]));
	}
	return j;
}

double portcullis::Junction::calcCodingPotential(GenomeMapper& gmap, KmerMarkovModel& exon, KmerMarkovModel& intron) {
	const char* ref = this->intron->ref.name.c_str();
	const bool neg = getConsensusStrand() == Strand::NEGATIVE;
	string left_exon = gmap.fetchBases(ref, this->intron->start - 82, this->intron->start - 2);
	if (neg) {
		left_exon = SeqUtils::reverseComplement(left_exon);
	}
	string left_intron = gmap.fetchBases(ref, this->intron->start, this->intron->start + 80);
	if (neg) {
		left_intron = SeqUtils::reverseComplement(left_intron);
	}
	string right_intron = gmap.fetchBases(ref, this->intron->end - 80, this->intron->end);
	if (neg) {
		right_intron = SeqUtils::reverseComplement(right_intron);
	}
	string right_exon = gmap.fetchBases(ref, this->intron->end + 1, this->intron->end + 81);
	if (neg) {
		right_exon = SeqUtils::reverseComplement(right_exon);
	}
	/*
	cout << "Left exon   : " << this->consensusStrand << " : " << left_exon << endl;
	cout << "Left intron : " << this->consensusStrand << " : " << left_intron << endl;
	cout << "Right intron: " << this->consensusStrand << " : " << right_intron << endl;
	cout << "Right exon  : " << this->consensusStrand << " : " << right_exon << endl;
	 */
	this->codingPotential = (exon.getScore(left_exon) - intron.getScore(left_exon))
							+ (intron.getScore(left_intron) - exon.getScore(left_intron))
							+ (intron.getScore(right_intron) - exon.getScore(right_intron))
							+ (exon.getScore(right_exon) - intron.getScore(right_exon));
	return this->codingPotential;
}

portcullis::SplicingScores portcullis::Junction::calcSplicingScores(GenomeMapper& gmap, KmerMarkovModel& donorT, KmerMarkovModel& donorF,
		KmerMarkovModel& acceptorT, KmerMarkovModel& acceptorF,
		PosMarkovModel& donorP, PosMarkovModel& acceptorP) {
	const char* ref = this->intron->ref.name.c_str();
	const bool neg = getConsensusStrand() == Strand::NEGATIVE;
	string left = gmap.fetchBases(ref, intron->start - 3, intron->start + 20);
	if (neg) {
		left = SeqUtils::reverseComplement(left);
	}
	string right = gmap.fetchBases(ref, intron->end - 20, intron->end + 2);
	if (neg) {
		right = SeqUtils::reverseComplement(right);
	}
	string donorseq = neg ? right : left;
	string acceptorseq = neg ? left : right;
	SplicingScores ss;
	ss.positionWeighting = donorP.getScore(donorseq) + acceptorP.getScore(acceptorseq);
	ss.splicingSignal = (donorT.getScore(donorseq) - donorF.getScore(donorseq))
						+ (acceptorT.getScore(acceptorseq) - acceptorF.getScore(acceptorseq));
	this->setPositionWeightScore(ss.positionWeighting);
	this->setSplicingSignal(ss.splicingSignal);
	return ss;
}

double portcullis::Junction::calcJunctionAnchorDepthLogDeviation(size_t i) const {
	double Ni = junctionAnchorDepth[i]; // Actual count at this position
	if (Ni == 0.0) Ni = 0.000000000001; // Ensure some value > 0 here otherwise we get -infinity later.
	double Pi = 1.0 - ((double) i / (double) (this->meanReadLength / 2.0)); // Likely scale at this position
	double Ei = (double) this->getNbSplicedAlignments() * Pi; // Expected count at this position
	double Xi = log2(Ni / Ei);
	return Xi;
}
