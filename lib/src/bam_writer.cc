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

#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::make_shared;
using std::shared_ptr;
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
#include <htslib/bgzf.h>

#include <portcullis/bam/bam_alignment.hpp>
using portcullis::bam::BamAlignment;

#include <portcullis/bam/bam_writer.hpp>

void portcullis::bam::BamWriter::open(bam_hdr_t* header) {
	// split
	fp = bgzf_open(bamFile.c_str(), "w");
	if (fp == NULL) {
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Could not open output BAM file: ") + bamFile.string()));
	}
	if (bam_hdr_write(fp, header) != 0) {
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Could not write header into: ") + bamFile.string()));
	}
}

int portcullis::bam::BamWriter::write(const BamAlignment& ba) {
	return bam_write1(fp, ba.getRaw());
}

void portcullis::bam::BamWriter::close() {
	bgzf_close(fp);
}

