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

#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using std::unique_ptr;
using std::string;
using std::unordered_map;
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
#include <htslib/bgzf.h>

#include <portcullis/bam/bam_alignment.hpp>
using portcullis::bam::BamAlignment;
using portcullis::bam::BamAlignmentPtr;
using portcullis::bam::RefSeqPtr;
using portcullis::bam::RefSeqPtrList;
using portcullis::bam::RefSeqPtrIndexMap;


namespace portcullis {
namespace bam {


class BamReader {

private:

	path bamFile;
	uint16_t threads;

	BGZF *fp;
	bam_hdr_t* header;
	bam1_t* c;
	hts_idx_t* index;
	hts_itr_t * iter;

	BamAlignment b;

public:

	BamReader(const path& _bamFile);

	virtual ~BamReader();

	// **** Methods for extracting reference target sequences ********

	shared_ptr<RefSeqPtrList> createRefList();

	shared_ptr<RefSeqPtrIndexMap> createRefMap();

	shared_ptr<RefSeqPtrIndexMap> createRefMap(const RefSeqPtrList& refList);

	static inline int bam_iter_read(BGZF* fp, hts_itr_t* iter, bam1_t *b) { return iter ? hts_itr_next(fp, iter, b, 0) : bam_read1(fp, b); }

	bam_hdr_t* getHeader() const { return header; }

	string bamDetails() const;

	void open();

	void close();

	bool next();

	const BamAlignment& current() const;

	void setRegion(const int32_t seqIndex, const int32_t start, const int32_t end);

	bool isCoordSortedBam();
};

}
}
