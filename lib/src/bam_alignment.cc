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
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::cout;
using std::endl;
using std::shared_ptr;
using std::make_shared;
using std::string;
using std::vector;
using std::stringstream;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::lexical_cast;

#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <portcullis/bam/bam_alignment.hpp>
using portcullis::bam::CigarOp;

portcullis::bam::CigarOp::CigarOp(const string& cigar) {
	if (cigar.size() < 2) {
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Can't extract cigar op from string when string length is < 2: ")
							  + cigar));
	}
	type = cigar[cigar.size() - 1];
	string len = cigar.substr(0, cigar.size() - 1);
	length = lexical_cast<uint32_t>(len);
}

vector<CigarOp> portcullis::bam::CigarOp::createFullCigarFromString(const string& cigar) {
	vector<CigarOp> c;
	uint32_t i = 0;
	while (i < cigar.size()) {
		size_t next = cigar.find_first_not_of("0123456789", i + 1);
		string op;
		if (next != string::npos) {
			op = cigar.substr(i, next - i + 1);
			i = next + 1;
		}
		c.push_back(CigarOp(op));
	}
	return c;
}

void portcullis::bam::BamAlignment::init() {
	alFlag = b->core.flag;
	position = b->core.pos;
	refId = b->core.tid;
	mateId = b->core.mtid;
	matePos = b->core.mpos;
	uint32_t* c = bam_get_cigar(b);
	alignedLength = 0;
	// Record the cigar in an easy to access format and calculate the aligned length
	// along the way
	cigar.clear();
	for (uint32_t i = 0; i < b->core.n_cigar; i++) {
		CigarOp op(bam_cigar_opchr(c[i]), bam_cigar_oplen(c[i]));
		cigar.push_back(op);
		if (CigarOp::opConsumesReference(op.type)) {
			alignedLength += op.length;
		}
	}
	// Determine actual read strand
	// Try deriving from XS tag first if it's present.  If not then try to
	// work it all out from the protocol suggested and the reverse strand
	// flag.
	Strand s = getXSStrand();
	if (s != Strand::UNKNOWN) {
		strand = s;
	}
	else {
		strand = calcStrand();
	}
}

portcullis::bam::Strand portcullis::bam::BamAlignment::calcStrand() {
	Strand strand = Strand::UNKNOWN;
	if (strandedness == Strandedness::FIRSTSTRAND) {
		if (orientation == Orientation::FR) {
			// Second mate should have correct strand (flip the strand for the first mate)
			if (isFirstMate()) {
				strand = isReverseStrand() ? Strand::POSITIVE : Strand::NEGATIVE;
			}
			else {
				strand = isReverseStrand() ? Strand::NEGATIVE : Strand::POSITIVE;
			}
		}
		else if (orientation == Orientation::RF) {
			if (isFirstMate()) {
				strand = isReverseStrand() ? Strand::NEGATIVE : Strand::POSITIVE;
			}
			else {
				strand = isReverseStrand() ? Strand::POSITIVE : Strand::NEGATIVE;
			}
		}
		else if (orientation == Orientation::F || orientation == Orientation::FF) {
			strand = isReverseStrand() ? Strand::POSITIVE : Strand::NEGATIVE;
		}
		else if (orientation == Orientation::R || orientation == Orientation::RR) {
			strand = isReverseStrand() ? Strand::NEGATIVE : Strand::POSITIVE;
		}
	}
	else if (strandedness == Strandedness::SECONDSTRAND) {
		if (orientation == Orientation::FR) {
			// First mate should have correct strand (flip the strand for the second mate)
			if (isFirstMate()) {
				strand = isReverseStrand() ? Strand::NEGATIVE : Strand::POSITIVE;
			}
			else {
				strand = isReverseStrand() ? Strand::POSITIVE : Strand::NEGATIVE;
			}
		}
		else if (orientation == Orientation::RF) {
			if (isFirstMate()) {
				strand = isReverseStrand() ? Strand::POSITIVE : Strand::NEGATIVE;
			}
			else {
				strand = isReverseStrand() ? Strand::NEGATIVE : Strand::POSITIVE;
			}
		}
		else if (orientation == Orientation::F || orientation == Orientation::FF) {
			strand = isReverseStrand() ? Strand::NEGATIVE : Strand::POSITIVE;
		}
		else if (orientation == Orientation::R || orientation == Orientation::RR) {
			strand = isReverseStrand() ? Strand::POSITIVE : Strand::NEGATIVE;
		}
	}
	return strand;
}

/**
 * Creates an empty samtools bam alignment
 */
portcullis::bam::BamAlignment::BamAlignment() {
	b = bam_init1();
	managed = true;
	alFlag = 0;
	position = -1;
	matePos = -1;
	alignedLength = 0;
	refId = -1;
	mateId = -1;
	strandedness = Strandedness::UNKNOWN;
	orientation = Orientation::UNKNOWN;
}

/**
 * Makes a deep copy of the provided samtools bam alignment
 * @param _b Samtools bam alignment
 */
portcullis::bam::BamAlignment::BamAlignment(bam1_t* _b, bool duplicate, Strandedness _strandedness, Orientation _orientation) {
	b = duplicate ? bam_dup1(_b) : _b;
	managed = duplicate;
	strandedness = _strandedness;
	orientation = _orientation;
	init();
}

/**
 * Makes a deep copy of an existing BamAlignment
 * @param other
 */
portcullis::bam::BamAlignment::BamAlignment(const shared_ptr<BamAlignment> other) {
	b = bam_dup1(other->b);
	managed = true;
	strandedness = other->strandedness;
	orientation = other->orientation;
	init();
}

/**
 * Makes a deep copy of an existing BamAlignment
 * @param other
 */
portcullis::bam::BamAlignment::BamAlignment(const BamAlignment& other) {
	b = bam_dup1(other.b);
	managed = true;
	strandedness = other.strandedness;
	orientation = other.orientation;
	init();
}

/**
 * Deletes the underlying samtools bam alignment only if it is managed (owned)
 * by this object
 */
portcullis::bam::BamAlignment::~BamAlignment() {
	if (managed)
		bam_destroy1(b);
}

void portcullis::bam::BamAlignment::setRaw(bam1_t* b) {
	this->b = b;
	managed = false;
	init();
}

bam1_t* portcullis::bam::BamAlignment::getRaw() const {
	return b;
}

/**
 * Looks for an XS tag in the alignment, returns the strand if found, or '?'
 * @return
 */
portcullis::bam::Strand portcullis::bam::BamAlignment::getXSStrand() const {
	char xs[2] = {'X', 'S'};
	uint8_t* res = bam_aux_get(b, xs);
	char c = res != 0 ? bam_aux2A(res) : '?';
	return strandFromChar(c);
}

string portcullis::bam::BamAlignment::deriveName() const {
	string qName = bam_get_qname(b);
	return isPaired() ?
		   qName + (this->isFirstMate() ?
					"_R1" :
					this->isSecondMate() ?
					"_R2" :
					"_R?") :
			   qName;
}

string portcullis::bam::BamAlignment::getQuerySeq() const {
	stringstream ss;
	for (int32_t i = 0; i < b->core.l_qseq; ++i) {
		ss << seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
	}
	return ss.str();
}

string portcullis::bam::BamAlignment::getQuerySeqAfterClipping() const {
	return getQuerySeqAfterClipping(getQuerySeq());
}

string portcullis::bam::BamAlignment::getQuerySeqAfterClipping(const string& seq) const {
	int32_t start = getStart();
	int32_t end = getEnd();
	int32_t clippedStart = cigar.front().type == BAM_CIGAR_SOFTCLIP_CHAR ? start + cigar.front().length : start;
	int32_t clippedEnd = cigar.back().type == BAM_CIGAR_SOFTCLIP_CHAR ? end - cigar.back().length : end;
	int32_t deltaStart = clippedStart - start;
	int32_t deltaEnd = end - clippedEnd;
	return seq.substr(deltaStart, seq.size() - deltaStart - deltaEnd + 1);
}

/**
 * We assume that orientatio
 * @param orientation
 * @return
 */
bool portcullis::bam::BamAlignment::calcIfProperPair(Orientation orientation) const {
	if (!this->isPaired() || !this->isMateMapped()) {
		return false;
	}
	if (this->refId != this->mateId) {
		return false;
	}
	bool diffStrand = this->isReverseStrand() != this->isMateReverseStrand();
	bool posGap = !this->isReverseStrand() ? this->position < this->matePos : this->position > this->matePos;
	if (orientation == Orientation::FR) {
		return diffStrand && posGap;
	}
	else if (orientation == Orientation::RF) {
		return diffStrand && !posGap;
	}
	else if (orientation == Orientation::FF) {
		return !diffStrand && posGap;
	}
	else if (orientation == Orientation::RR) {
		return !diffStrand && !posGap;
	}
	else {
		return false;
	}
}

bool portcullis::bam::BamAlignment::isSplicedRead() const {
	for (const auto & op : cigar) {
		if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
			return true;
		}
	}
	return false;
}

uint32_t portcullis::bam::BamAlignment::getNbJunctionsInRead() const {
	int32_t nbJunctions = 0;
	for (const auto & op : cigar) {
		if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
			nbJunctions++;
		}
	}
	return nbJunctions;
}

uint32_t portcullis::bam::BamAlignment::calcNbAlignedBases(int32_t start, int32_t end, bool includeSoftClips) const {
	if (start > getEnd() || end < position) {
		string align = this->toString();
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Found an alignment that does not have a presence in the requested region.  Requested region: ") +
							  lexical_cast<string>(start) + "-" + lexical_cast<string>(end) + ".  Alignment: " + align));
	}
	int32_t count = 0;
	int32_t pos = position;
	for (const auto & op : cigar) {
		if (pos > end) {
			break;
		}
		// Include softclips in this calculation
		if (CigarOp::opConsumesReference(op.type) || (includeSoftClips && op.type == BAM_CIGAR_SOFTCLIP_CHAR)) {
			if (pos >= start) {
				count += op.length;
			}
			pos += op.length;
		}
	}
	return count;
}

string portcullis::bam::BamAlignment::getPaddedQuerySeq(int32_t start, int32_t end, int32_t& actual_start, int32_t& actual_end, const bool include_soft_clips) const {
	return getPaddedQuerySeq(this->getQuerySeq(), start, end, actual_start, actual_end, include_soft_clips);
}

string portcullis::bam::BamAlignment::getPaddedQuerySeq(const string& query_seq, int32_t start, int32_t end, int32_t& actual_start, int32_t& actual_end, const bool include_soft_clips) const {
	if (start > getEnd() || end < position)
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Found an alignment that does not have a presence in the requested region")));
	int32_t qPos = 0;
	int32_t rPos = position;
	string query = include_soft_clips ? query_seq : this->getQuerySeqAfterClipping(query_seq);
	stringstream ss;
	for (const auto & op : cigar) {
		bool consumesRef = CigarOp::opConsumesReference(op.type);
		bool consumesQuery = CigarOp::opConsumesQuery(op.type) && (include_soft_clips || op.type != BAM_CIGAR_SOFTCLIP_CHAR);
		// Skips any cigar ops before start position
		if (rPos < start) {
			if (consumesRef) rPos += op.length;
			if (consumesQuery) qPos += op.length;
			continue;
		}
		// Stop once we get to the end of the region, and make sure we don't end on a refskip that exceeds our limit
		if (((rPos > end && op.type != BAM_CIGAR_INS_CHAR) || (op.type == BAM_CIGAR_REFSKIP_CHAR && rPos + op.length > end))) break;
		if (consumesQuery) {
			// Don't return anything that runs off the end cap
			int32_t len = rPos + op.length > end && op.type != BAM_CIGAR_INS_CHAR ? end - rPos + 1 : op.length;
			if (len == 0) {
				BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
										  "Can't extract cigar op sequence from query string when length has been calculated as 0.")
									  + "\nLimits: " + lexical_cast<string>(start) + "," + lexical_cast<string>(end)
									  + "\nQuery sequence: " + query + " (" + lexical_cast<string>(query.size()) + ")"
									  + "\nCigar: " + getCigarAsString()
									  + "\nCurrent op: " + op.toString()
									  + "\nCurrent position in query: " + lexical_cast<string>(qPos)
									  + "\nRequested length: " + lexical_cast<string>(len)
									  + "\nCurrent position in reference: " + lexical_cast<string>(rPos)
									  + "\nCurrent output: " + ss.str()));
			}
			//cout << qPos << endl;
			if (qPos < 0 || qPos + len > (int32_t)query.size()) {
				BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
										  "Can't extract cigar op sequence from query string.")
									  + "\nLimits: " + lexical_cast<string>(start) + "," + lexical_cast<string>(end)
									  + "\nQuery sequence: " + query + " (" + lexical_cast<string>(query.size()) + ")"
									  + "\nCigar: " + getCigarAsString()
									  + "\nCurrent op: " + op.toString()
									  + "\nCurrent position in query: " + lexical_cast<string>(qPos)
									  + "\nRequested length: " + lexical_cast<string>(len)
									  + "\nCurrent position in reference: " + lexical_cast<string>(rPos)
									  + "\nCurrent output: " + ss.str()));
			}
			ss << query.substr(qPos, len);
		}
		else if (consumesRef) {   // i.e. consumes reference but not query (DEL or REF_SKIP ops)
			uint32_t len = rPos + op.length > end ? end - rPos + 1 : op.length; // Make sure we don't exceed our ref end limit
			string s;
			s.resize(len);
			std::fill(s.begin(), s.end(), BAM_CIGAR_DIFF_CHAR);
			ss << s;
		}
		if (consumesRef) rPos += op.length;
		if (consumesQuery) qPos += op.length;
	}
	actual_start = position > start ? position : start;
	actual_end = rPos <= end ? rPos - 1 : end;
	return ss.str();
}

string portcullis::bam::BamAlignment::getPaddedGenomeSeq(const string& genomeSeq, int32_t start, int32_t end, int32_t q_start, int32_t q_end, const bool include_soft_clips) const {
	if (start > getEnd() || end < position)
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Found an alignment that does not have a presence in the requested region")));
	//int32_t pos = 0;
	//int32_t qPos = 0;
	int32_t rPos = position;
	int32_t startDelta = q_start - start;
	int32_t endDelta = end - q_end;
	if (startDelta < 0) {
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Query start position was before genomic region start position.  Query start: ") + lexical_cast<string>(q_start) + "; Genomic start: " + lexical_cast<string>(start)));
	}
	if (endDelta < 0) {
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Query end position was beyond genomic region end position.  Query end: ") + lexical_cast<string>(q_end) + "; Genomic end: " + lexical_cast<string>(end)));
	}
	stringstream ss;
	for (const auto & op : cigar) {
		bool consumesRef = CigarOp::opConsumesReference(op.type);
		bool consumesQuery = CigarOp::opConsumesQuery(op.type) && (include_soft_clips || op.type != BAM_CIGAR_SOFTCLIP_CHAR);
		// Skips any cigar ops before start position
		if (rPos < q_start) {
			if (consumesRef) rPos += op.length;
			//if (consumesQuery) qPos += op.length;
			continue;
		}
		// Ends cigar loop once we reach the end of the region (unless we have an insertion op here... then proceed)
		if (rPos > q_end && op.type != BAM_CIGAR_INS_CHAR) break;
		if (consumesRef) {
			int32_t seqOffset = rPos - start;
			int32_t len = rPos + op.length > q_end ? q_end - rPos + 1 : op.length;
			if (seqOffset < 0 || seqOffset + len > (int32_t)genomeSeq.size()) {
				BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
										  "Can't extract cigar op sequence from extracted genome region.\nCurrent position in extracted genome region: ")
									  + lexical_cast<string>(seqOffset) +
									  "; Genome region length: " + lexical_cast<string>(genomeSeq.size()) +
									  "\nFull cigar: " + this->getCigarAsString() +
									  "\nCurrent cigar op: " + op.type + lexical_cast<string>(op.length) +
									  "\nCurrent genomic position: " + lexical_cast<string>(rPos) +
									  "\nAlignment region: " + lexical_cast<string>(position) + "," + lexical_cast<string>(position + alignedLength) +
									  "\nGenome region: " + lexical_cast<string>(start) + "," + lexical_cast<string>(end) +
									  "\nQuery region: " + lexical_cast<string>(q_start) + "," + lexical_cast<string>(q_end)));
			}
			//cout << seqOffset << endl;
			ss << genomeSeq.substr(seqOffset, len);
		}
		else if (consumesQuery) {   // Consumes query but not reference ('I' op)
			string s;
			s.resize(op.length);
			std::fill(s.begin(), s.end(), BAM_CIGAR_DIFF_CHAR);
			ss << s;
		}
		if (consumesRef) rPos += op.length;
		//if (consumesQuery) qPos += op.length;
	}
	return ss.str();
}

string portcullis::bam::BamAlignment::toString() const {
	return toString(false);
}

string portcullis::bam::BamAlignment::toString(bool afterClipping) const {
	uint32_t start = afterClipping && cigar.front().type == BAM_CIGAR_SOFTCLIP_CHAR ? position + cigar.front().length : position;
	uint32_t end = afterClipping && cigar.back().type == BAM_CIGAR_SOFTCLIP_CHAR ? getEnd() - cigar.back().length : getEnd();
	stringstream ss;
	ss << refId << "(" << start << "-" << end << ")" << (this->isReverseStrand() ? "-" : "+");
	return ss.str();
}