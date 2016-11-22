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

#include <ctime>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::time_t;
using std::difftime;
using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;
using std::stringstream;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::filesystem::last_write_time;
using boost::lexical_cast;

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <portcullis/bam/bam_master.hpp>


bool portcullis::bam::BamHelper::isCoordSortedBam(const path& bamFile) {
	BGZF *fp;
	// split
	fp = strcmp(bamFile.c_str(), "-") ? bgzf_open(bamFile.c_str(), "r") : bgzf_dopen(fileno(stdin), "r");
	if (fp == NULL) {
		BOOST_THROW_EXCEPTION(BamException() << BamErrorInfo(string(
								  "Could not open input BAM files: ") + bamFile.string()));
	}
	bam_hdr_t* header = bam_hdr_read(fp);
	string headerText = header->text;
	bool found = false;
	if (headerText.find("SO:coordinate") != std::string::npos) {
		found = true;
	}
	bam_hdr_destroy(header);
	return found;
}

bool portcullis::bam::BamHelper::isNewerIndexPresent(const path& bamFile, bool useCsi) {
	path idxFile = path(bamFile.string() + (useCsi ? ".csi" : ".bai"));
	if (!exists(idxFile))
		return false;
	time_t bamTime = last_write_time( bamFile ) ;
	time_t idxTime = last_write_time( bamFile ) ;
	return difftime(bamTime, idxTime) > 0;
}

/**
 * Creates a command that can be used to merge multiple BAM files with samtools
 * @param samtoolsExe The path to samtools
 * @param bamFiles The paths to each BAM file to merge
 * @param mergedBamFile The output file
 * @param threads Number of threads to use during merging
 * @return Command line
 */
string portcullis::bam::BamHelper::createMergeBamCmd(const vector<path>& bamFiles,
		const path& mergedBamFile,
		uint16_t threads) {
	stringstream inputFiles;
	for (path p : bamFiles) {
		inputFiles << " " << p.c_str();
	}
	return string("samtools merge -f -@ ") + lexical_cast<string>(threads) +
		   " " + mergedBamFile.string() +
		   inputFiles.str();
}


/**
 * Creates a samtools command that can be used to sort a bam file
 * @param samtoolsExe The path to samtools
 * @param unsortedFile The bam file that needs sorting
 * @param sortedFile The path to the new sorted bam file which will be created
 * @param sortByName If true, bam entries are sorted by name, otherwise by position
 * @param threads Number of threads to use
 * @param memory Amount of memory to request
 * @return The command that can be used to sort the bam file
 */
string portcullis::bam::BamHelper::createSortBamCmd(const path& unsortedFile,
		const path& sortedFile,
		bool sortByName,
		uint16_t threads,
		const string& memory) {
	return string("samtools sort -@ ") + lexical_cast<string>(threads) +
		   " -m " + memory + " " + (sortByName ? "-n " : "") + unsortedFile.string() +
		   " " + sortedFile.string();
}

/**
 * Creates a samtools command that can be used to index a sorted bam file
 * @param samtoolsExe The path to samtools
 * @param sortedBam Path to a sorted bam file to index
 * @return The command that can be used to index the sorted bam file
 */
string portcullis::bam::BamHelper::createIndexBamCmd(const path& sortedBam, bool useCsi) {
	return string("samtools index ") + (useCsi ? "-c " : "") + sortedBam.string();
}



