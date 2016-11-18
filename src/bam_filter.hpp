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

#include <portcullis/bam/bam_alignment.hpp>
using portcullis::bam::BamAlignment;
using portcullis::bam::BamAlignmentPtr;

#include <portcullis/junction_system.hpp>


namespace portcullis {
    
typedef boost::error_info<struct BamFilterError,string> BamFilterErrorInfo;
struct BamFilterException: virtual boost::exception, virtual std::exception { };

enum class ClipMode {
    HARD,
    SOFT,
    COMPLETE
};

static string clipToString(ClipMode cm) {
    
    switch(cm) {
        case ClipMode::HARD:
            return "HARD";
        case ClipMode::SOFT:
            return "SOFT";
        case ClipMode::COMPLETE:
            return "COMPLETE";
    }

    return "COMPLETE";
}

static ClipMode clipFromString(string cm) {
    
    if (boost::iequals(cm, "HARD")) {
        return ClipMode::HARD;
    }
    else if (boost::iequals(cm, "SOFT")) {
        return ClipMode::SOFT;
    }
    else if (boost::iequals(cm, "COMPLETE")) {
        return ClipMode::COMPLETE;
    }
    
    BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Unrecognised clip mode: ") + cm));
}

class BamFilter {

private:
    
    path junctionFile;
    path bamFile;
    path outputBam;
    Strandedness strandSpecific;
    Orientation orientation;
    ClipMode clipMode;
    bool saveMSRs;
    bool useCsi;
    bool verbose;
    
public:
    
    BamFilter(const path& _junctionFile, const path& _bamFile, const path& _outputBam);
    
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
    bool containsJunctionInSystem(const BamAlignment& al, const RefSeqPtrList& refs, JunctionSystem& js);
    
    BamAlignmentPtr clipMSR(const BamAlignment& al, const RefSeqPtrList& refs, JunctionSystem& js, bool& allBad);
       

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

    Strandedness getStrandSpecific() const {
        return strandSpecific;
    }

    void setStrandSpecific(Strandedness strandSpecific) {
        this->strandSpecific = strandSpecific;
    }
    
    Orientation getOrientation() const {
        return orientation;
    }

    void setOrientation(Orientation orientation) {
        this->orientation = orientation;
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
    
    bool isUseCsi() const {
        return useCsi;
    }

    void setUseCsi(bool useCsi) {
        this->useCsi = useCsi;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }




    
    void filter();
  
    static string title() {
        return string("Portcullis BAM Filter Mode Help.");
    }
    
    static string description() {
        return string("Removes alignments associated with bad junctions from BAM file");
    }
    
    static string usage() {
        return string("portcullis bamfilt [options] <junction-file> <bam-file>");
    }
    
    
    static int main(int argc, char *argv[]);
};
}