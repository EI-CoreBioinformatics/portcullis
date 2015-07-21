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

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
namespace bfs = boost::filesystem;
namespace po = boost::program_options;

#include "bam/bam_alignment.hpp"
using portcullis::bam::BamAlignment;
using portcullis::bam::BamAlignmentPtr;

#include "junction_system.hpp"


namespace portcullis {
    
typedef boost::error_info<struct BamFilterError,string> BamFilterErrorInfo;
struct BamFilterException: virtual boost::exception, virtual std::exception { };

enum ClipMode {
    HARD,
    SOFT,
    COMPLETE
};

static string clipToString(ClipMode cm) {
    
    switch(cm) {
        case HARD:
            return "HARD";
        case SOFT:
            return "SOFT";
        case COMPLETE:
            return "COMPLETE";
    }

    return "COMPLETE";
}

static ClipMode clipFromString(string cm) {
    
    if (boost::iequals(cm, "HARD")) {
        return HARD;
    }
    else if (boost::iequals(cm, "SOFT")) {
        return SOFT;
    }
    else if (boost::iequals(cm, "COMPLETE")) {
        return COMPLETE;
    }
    
    BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Unrecognised clip mode: ") + cm));
}

class BamFilter {

private:
    
    path junctionFile;
    path bamFile;
    path outputBam;
    StrandSpecific strandSpecific;
    ClipMode clipMode;
    bool saveMSRs;
    bool verbose;
    
public:
    
    BamFilter(const path& _junctionFile, const path& _bamFile, const path& _outputBam, bool _verbose);
    
    virtual ~BamFilter() {
    }
    
    
protected:
    
    /**
     * Checks a given alignment to see if it exists in the given junction system
     * @param al Alignment to check
     * @param refs References
     * @param js The junction system containing good junctions to keep
     * @return Whether or not the alignment contains a junction found in the junction system
     */
    bool containsJunctionInSystem(const BamAlignment& al, vector<RefSeq>& refs, JunctionSystem& js);
    
    BamAlignmentPtr clipMSR(const BamAlignment& al, vector<RefSeq>& refs, JunctionSystem& js, bool& allBad);
       

public:
    
    path getBamFile() const {
        return bamFile;
    }

    void setBamFile(path bamFile) {
        this->bamFile = bamFile;
    }

    path getJunctionFile() const {
        return junctionFile;
    }

    void setJunctionFile(path junctionFile) {
        this->junctionFile = junctionFile;
    }

    path getOutputBam() const {
        return outputBam;
    }

    void setOutputBam(path outputBam) {
        this->outputBam = outputBam;
    }

    StrandSpecific getStrandSpecific() const {
        return strandSpecific;
    }

    void setStrandSpecific(StrandSpecific strandSpecific) {
        this->strandSpecific = strandSpecific;
    }
    
    ClipMode getClipMode() const {
        return clipMode;
    }

    void setClipMode(ClipMode clipMode) {
        this->clipMode = clipMode;
    }
    
    bool isSaveMSRs() const {
        return saveMSRs;
    }

    void setSaveMSRs(bool saveMSRs) {
        this->saveMSRs = saveMSRs;
    }



    
    void filter();
  
    static string helpMessage() {
        return string("\nPortcullis BAM Filter Mode Help.\n\n") +
                      "Removes alignments associated with bad junctions from BAM file\n\n" + 
                      "Usage: portcullis bamfilt [options] <junction-file> <bam-file>\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]);
};
}