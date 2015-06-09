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
#include <iostream>
#include <vector>
#include <memory>
using std::boolalpha;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using namespace boost::filesystem;
using boost::filesystem::path;
namespace po = boost::program_options;

#include "samtools_helper.hpp"
#include "intron.hpp"
#include "junction.hpp"
#include "junction_system.hpp"
#include "prepare.hpp"
using portcullis::BamAlignment;
using portcullis::BamReader;
using portcullis::GenomeMapper;
using portcullis::Intron;
using portcullis::Junction;
using portcullis::JunctionSystem;



namespace portcullis {

const string DEFAULT_JUNC_OUTPUT_DIR = "portcullis_junc_out";
const string DEFAULT_JUNC_OUTPUT_PREFIX = "portcullis";
const uint16_t DEFAULT_JUNC_THREADS = 1;

typedef boost::error_info<struct JunctionBuilderError,string> JunctionBuilderErrorInfo;
struct JunctionBuilderException: virtual boost::exception, virtual std::exception { };

class JunctionBuilder {
private:

    // Can set these from the outside via the constructor
    PreparedFiles* prepData;
    path outputDir;
    string outputPrefix;
    StrandSpecific strandSpecific;
    uint16_t threads;
    bool fast;
    path samtoolsExe;
    bool verbose;
    
    // The set of distinct junctions found in the BAM file
    JunctionSystem junctionSystem;
    SplicedAlignmentMap splicedAlignmentMap;
    
    

protected:
    
    path getUnsplicedBamFile() {
        return path(outputDir.string() + "/" + outputPrefix + ".unspliced.bam");
    }
    
    path getSplicedBamFile() {
        return path(outputDir.string() + "/" + outputPrefix + ".spliced.bam");
    }
    
    path getAssociatedIndexFile(path bamFile) {
        return path(bamFile.string() + ".bai");
    }
        

public:

    
    JunctionBuilder(const path& _prepDir, const path& _outputDir, string _outputPrefix, uint16_t _threads, bool _fast, bool _verbose);
    
    virtual ~JunctionBuilder();
    
    
    
    /**
     * Populates the set of distinct junctions.  
     * 
     * Also outputs all the unspliced alignments to a separate file if requested
     */
    void process();
    
    path getSamtoolsExe() const {
        return samtoolsExe;
    }

    void setSamtoolsExe(path samtoolsExe) {
        this->samtoolsExe = samtoolsExe;
    }

    
    static string helpMessage() {
        return string("\nPortcullis Junction Builder Mode Help.\n\n") +
                      "Usage: portcullis junc [options] <prep_data_dir> \n\n" +
                      "Run \"portcullis prep ...\" to generate data suitable for junction finding before running \"portcullis junc ...\"\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]);
};
}

