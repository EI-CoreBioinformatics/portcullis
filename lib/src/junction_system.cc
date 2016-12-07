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

#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::shared_ptr;
using std::unordered_set;

#include <boost/lexical_cast.hpp>
#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
using boost::lexical_cast;
using boost::timer::auto_cpu_timer;

#include <portcullis/bam/depth_parser.hpp>
using portcullis::bam::DepthParser;

#include <portcullis/intron.hpp>
#include <portcullis/junction.hpp>
#include <portcullis/seq_utils.hpp>
using portcullis::Intron;
using portcullis::IntronHasher;
using portcullis::Junction;
using portcullis::JunctionPtr;
using portcullis::SeqUtils;

#include <portcullis/junction_system.hpp>

size_t portcullis::JunctionSystem::createJunctionGroup(size_t index, vector<JunctionPtr>& group) {
	JunctionPtr junc = junctionList[index];
	group.push_back(junc);
	bool foundMore = false;
	for (size_t j = index + 1; j < junctionList.size(); j++) {
		JunctionPtr next = junctionList[j];
		if (junc->sharesDonorOrAcceptor(next)) {
			group.push_back(next);
			junc = next;
		}
		else {
			return j - 1;
		}
	}
	return foundMore ? index : junctionList.size() - 1;
}

void portcullis::JunctionSystem::findJunctions(const int32_t refId, JunctionList& subset) {
	subset.clear();
	for (JunctionPtr j : junctionList) {
		if (j->getIntron()->ref.index == refId) {
			subset.push_back(j);
		}
	}
}

portcullis::JunctionSystem::JunctionSystem() {
	minQueryLength = 0;
	meanQueryLength = 0.0;
	maxQueryLength = 0;
	distinctJunctions.clear();
	junctionList.clear();
}

portcullis::JunctionSystem::JunctionSystem(shared_ptr<vector<RefSeqPtr>> refs) : JunctionSystem() {
	this->refs = refs;
}

portcullis::JunctionSystem::JunctionSystem(path junctionFile) : JunctionSystem() {
	load(junctionFile);
}

portcullis::JunctionSystem::JunctionSystem(JunctionList& jl) : JunctionSystem() {
	for (auto & j : jl) {
		this->addJunction(j);
	}
}

portcullis::JunctionSystem::~JunctionSystem() {
	distinctJunctions.clear();
	junctionList.clear();
}

const JunctionList& portcullis::JunctionSystem::getJunctions() const {
	return junctionList;
}

size_t portcullis::JunctionSystem::size() {
	//assert(distinctJunctions.size() == junctionList.size());
	return distinctJunctions.size();
}

void portcullis::JunctionSystem::setQueryLengthStats(int32_t min, double mean, int32_t max) {
	this->minQueryLength = min;
	this->meanQueryLength = mean;
	this->maxQueryLength = max;
}

void portcullis::JunctionSystem::addJunction(JunctionPtr j) {
	j->clearAlignments();
	distinctJunctions[*(j->getIntron())] = j;
	junctionList.push_back(j);
}

/**
 * Appends a new copy of all the junctions in the other junction system to this
 * junction system, without the Bam alignments associated with them.
 * @param other The other junctions system containing junctions to be added to this.
 */
void portcullis::JunctionSystem::append(JunctionSystem& other) {
	for (const auto & j : other.getJunctions()) {
		this->addJunction(j);
	}
}

bool portcullis::JunctionSystem::addJunctions(const BamAlignment& al, const size_t startOp, const int32_t offset) {
	bool foundJunction = false;
	const size_t nbOps = al.getNbCigarOps();
	const int32_t refId = al.getReferenceId();
	int32_t lStart = offset;
	int32_t lEndExc = lStart; // End of left anchor exclusive (i.e. +1)
	int32_t rStart = lStart;
	int32_t rEndExc = lStart; // End of right anchor exclusive (i.e. +1)
	for (size_t i = startOp; i < nbOps; i++) {
		CigarOp op = al.getCigarOpAt(i);
		if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
			foundJunction = true;
			const int32_t refLength = refs->at(refId)->length;
			rStart = lEndExc + op.length;
			rEndExc = rStart;
			// Establish end position of right anchor in genomic coordinates
			size_t j = i + 1;
			while (j < nbOps
					&& rEndExc <= refLength
					&& al.getCigarOpAt(j).type != BAM_CIGAR_REFSKIP_CHAR) {
				CigarOp rOp = al.getCigarOpAt(j++);
				if (CigarOp::opConsumesReference(rOp.type)) {
					rEndExc += rOp.length;
				}
			}
			// Do some sanity checking... make sure there are not any strange N cigar ops that
			// drift over the edge of a reference sequence... seems like this can actually
			// happen in GSNAP!  I guess this is referring to reads that map over
			// target sequence boundaries
			if (rStart - 1 >= refLength) {
				rStart = refLength - 1;
			}
			if (rEndExc - 1 >= refLength) {
				rEndExc = refLength;
			}
			// Create the intron
			shared_ptr<Intron> location = make_shared<Intron>(
											  RefSeq(refId, refs->at(refId)->name, refLength),
											  lEndExc,
											  rStart - 1);
			// We should now have the complete junction location information
			JunctionMapIterator it = distinctJunctions.find(*location);
			// If we couldn't find this location in the hashmap, add a new
			// location / junction pair.  If we've seen this location before
			// then add this alignment to the existing junction
			if (it == distinctJunctions.end()) {
				JunctionPtr junction = make_shared<Junction>(location, lStart, rEndExc - 1);
				junction->addJunctionAlignment(al);
				distinctJunctions[*location] = junction;
				junctionList.push_back(junction);
			}
			else {
				JunctionPtr junction = it->second;
				junction->addJunctionAlignment(al);
				junction->extendAnchors(lStart, rEndExc - 1);
			}
			// Check if we have fully processed the cigar or not.  If not, then
			// that means that this cigar contains additional junctions, so
			// process those using recursion
			if (j < nbOps) {
				addJunctions(al, i + 1, rStart);
				break;
			}
		}
		else if (CigarOp::opConsumesReference(op.type)) {
			lEndExc += op.length;
		}
		// Ignore any other op types not already covered
	}
	return foundJunction;
}

void portcullis::JunctionSystem::findFlankingAlignments(const path& alignmentsFile) {
	auto_cpu_timer timer(1, " done. Wall time taken: %ws\n");
	// Maybe try to multi-thread this part
	BamReader reader(alignmentsFile);
	// Open the file
	reader.open();
	// Read the alignments around every junction and set appropriate metrics
	//size_t count = 0;
	//cout << endl;
	for (JunctionPtr j : junctionList) {
		j->processJunctionVicinity(
			reader,
			j->getIntron()->ref.length,
			maxQueryLength);
		//cout << count++ << endl;
	}
	reader.close();
}

void portcullis::JunctionSystem::calcCoverage(const path& alignmentsFile, Strandedness strandSpecific) {
	auto_cpu_timer timer(1, " done. Wall time taken: %ws\n");
	DepthParser dp(alignmentsFile, static_cast<uint8_t> (strandSpecific), false);
	vector<uint32_t> batch;
	while (dp.loadNextBatch(batch)) {
		JunctionList subset;
		findJunctions(dp.getCurrentRefIndex(), subset);
		for (JunctionPtr j : subset) {
			j->calcCoverage(batch);
		}
	}
}

void portcullis::JunctionSystem::calcMultipleMappingStats(SplicedAlignmentMap& map) {
	for (JunctionPtr j : junctionList) {
		j->calcMultipleMappingScore(map);
	}
}

void portcullis::JunctionSystem::calcJunctionStats() {
	if (junctionList.empty()) {
		return;
	}
	for (size_t i = 0; i < junctionList.size(); i++) {
		vector<JunctionPtr > junctionGroup;
		i = createJunctionGroup(i, junctionGroup);
		uint32_t maxReads = 0;
		size_t maxIndex = 0;
		bool uniqueJunction = junctionGroup.size() == 1;
		for (size_t j = 0; j < junctionGroup.size(); j++) {
			JunctionPtr junc = junctionGroup[j];
			if (maxReads < junc->getNbSplicedAlignments()) {
				maxReads = junc->getNbSplicedAlignments();
				maxIndex = j;
			}
			junc->setUniqueJunction(uniqueJunction);
		}
		junctionGroup[maxIndex]->setPrimaryJunction(true);
	}
	size_t i = 0;
	bool lastdiffseq = false;
	while (i < junctionList.size() - 1) {
		JunctionPtr first = junctionList[i];
		JunctionPtr second = junctionList[i + 1];
		int32_t diff = second->getIntron()->start - first->getIntron()->end;
		diff = diff < 0 ? 0 : diff;
		if (first->getIntron()->ref.index != second->getIntron()->ref.index) {
			first->setDistanceToNextUpstreamJunction(-1);
			second->setDistanceToNextDownstreamJunction(-1);
			if (i == 0 || lastdiffseq) {
				first->setDistanceToNextDownstreamJunction(-1);
			}
			if (i == junctionList.size() - 2) {
				second->setDistanceToNextUpstreamJunction(-1);
			}
			lastdiffseq = true;
		}
		else if (i == 0) {
			first->setDistanceToNextDownstreamJunction(-1);
			first->setDistanceToNextUpstreamJunction(diff);
			second->setDistanceToNextDownstreamJunction(diff);
			lastdiffseq = false;
		}
		else if (i == junctionList.size() - 2) {
			first->setDistanceToNextUpstreamJunction(diff);
			second->setDistanceToNextDownstreamJunction(diff);
			second->setDistanceToNextUpstreamJunction(-1);
			lastdiffseq = false;
		}
		else {
			first->setDistanceToNextUpstreamJunction(diff);
			second->setDistanceToNextDownstreamJunction(diff);
			lastdiffseq = false;
		}
		i++;
	}
	for (auto & junc : junctionList) {
		int32_t down = junc->getDistanceToNextDownstreamJunction();
		int32_t up = junc->getDistanceToNextUpstreamJunction();
		junc->setDistanceToNearestJunction(down == -1 || up == -1 ? max(down, up) : min(down, up));
		junc->setMeanReadLength(this->meanQueryLength);
		// Now we know the mean query length, confirm if this junction really is suspicious
		if (junc->isSuspicious()) {
			double prob = 1.0 - std::pow((junc->getMaxMMES() / (this->meanQueryLength / 2.0)), junc->getNbSplicedAlignments());
			if (prob > 0.99) {
				junc->setPotentialFalsePositive(true);
			}
		}
	}
}

void portcullis::JunctionSystem::sort() {
	std::sort(junctionList.begin(), junctionList.end(), JunctionComparator());
}

void portcullis::JunctionSystem::index() {
	for (size_t i = 0; i < this->size(); i++) {
		junctionList[i]->setId(i);
	}
}

void portcullis::JunctionSystem::saveAll(const path& outputPrefix, const string& source) {
	saveAll(outputPrefix, source, false, false, false);
}

void portcullis::JunctionSystem::saveAll(const path& outputPrefix, const string& source, bool bedscore, bool outputExonGFF, bool outputIntronGFF) {
	auto_cpu_timer timer(1, " = Wall time taken: %ws\n\n");
	string junctionReportPath = outputPrefix.string() + ".junctions.txt";
	string junctionFilePath = outputPrefix.string() + ".junctions.tab";
	string junctionGFFPath = outputPrefix.string() + ".junctions.exon.gff3";
	string intronGFFPath = outputPrefix.string() + ".junctions.intron.gff3";
	string junctionBEDAllPath = outputPrefix.string() + ".junctions.bed";
	/*cout << " - Saving junction report to: " << junctionReportPath << " ... ";
	cout.flush();

	// Print descriptive output to file
	ofstream junctionReportStream(junctionReportPath.c_str());
	outputDescription(junctionReportStream);
	junctionReportStream.close();

	cout << "done." << endl;*/
	cout << " - Saving junction table to: " << junctionFilePath << " ... ";
	cout.flush();
	// Print junction stats to file
	ofstream junctionFileStream(junctionFilePath.c_str());
	junctionFileStream << (*this) << endl;
	junctionFileStream.close();
	cout << "done." << endl;
	if (outputExonGFF) {
		cout << " - Saving junction GFF file to: " << junctionGFFPath << " ... ";
		cout.flush();
		// Print junction stats to file
		ofstream junctionGFFStream(junctionGFFPath.c_str());
		writeExonGFF(junctionGFFStream, source);
		junctionGFFStream.close();
		cout << "done." << endl;
	}
	if (outputIntronGFF) {
		cout << " - Saving intron GFF file to: " << intronGFFPath << " ... ";
		cout.flush();
		// Print junction stats to file
		ofstream intronGFFStream(intronGFFPath.c_str());
		writeIntronGFF(intronGFFStream, source);
		intronGFFStream.close();
		cout << "done." << endl;
	}
	// Output BED files
	cout << " - Saving BED file with all junctions to: " << junctionBEDAllPath << " ... ";
	cout.flush();
	// Print junctions in BED format to file
	outputBED(junctionBEDAllPath, CanonicalSS::ALL, source, bedscore);
	cout << "done." << endl;
}

void portcullis::JunctionSystem::outputDescription(std::ostream &strm) {
	for (JunctionPtr j : junctionList) {
		strm << "Junction " << j->getId() << ":" << endl;
		j->outputDescription(strm);
		strm << endl;
	}
}

void portcullis::JunctionSystem::writeExonGFF(std::ostream &strm, const string& source) {
	for (JunctionPtr j : junctionList) {
		j->outputJunctionGFF(strm, source);
	}
}

void portcullis::JunctionSystem::writeIntronGFF(std::ostream &strm, const string& source) {
	for (JunctionPtr j : junctionList) {
		j->outputIntronGFF(strm, source);
	}
}

void portcullis::JunctionSystem::outputBED(string& path, CanonicalSS type, const string& prefix, bool bedscore) {
	ofstream junctionBEDStream(path.c_str());
	outputBED(junctionBEDStream, type, prefix, bedscore);
	junctionBEDStream.close();
}

void portcullis::JunctionSystem::outputBED(std::ostream &strm, CanonicalSS type, const string& prefix, bool bedscore) {
	strm << "track name=\"junctions\" description=\"Portcullis V" << (version.empty() ? "X.X.X" : version) << " junctions\"" << endl;
	for (JunctionPtr j : junctionList) {
		if (type == CanonicalSS::ALL || j->getSpliceSiteType() == type) {
			j->outputBED(strm, prefix, bedscore);
		}
	}
}

void portcullis::JunctionSystem::load(const path& junctionTabFile) {
	load(junctionTabFile, false);
}

void portcullis::JunctionSystem::load(const path& junctionTabFile, const bool simple) {
	ifstream ifs(junctionTabFile.c_str());
	string line;
	// Loop through until end of file or we move onto the next ref seq
	while (std::getline(ifs, line)) {
		boost::trim(line);
		if (!line.empty() && line.find("index") == std::string::npos) {
			shared_ptr<Junction> j = Junction::parse(line);
			junctionList.push_back(j);
			if (!simple) {
				distinctJunctions[*(j->getIntron())] = j;
			}
		}
	}
	ifs.close();
}

JunctionPtr portcullis::JunctionSystem::getJunction(Intron& intron) const {
	try {
		return this->distinctJunctions.at(intron);
	}
	catch (std::out_of_range ex) {
		return nullptr;
	}
}

string portcullis::JunctionSystem::version = "";

