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
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::pair;
using std::map;
using std::unordered_map;
using std::unordered_set;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
using boost::timer::auto_cpu_timer;
using boost::filesystem::absolute;
using boost::filesystem::copy_file;
using boost::filesystem::remove;
using boost::filesystem::exists;
using boost::filesystem::create_symlink;
using boost::filesystem::create_directory;
using boost::filesystem::symbolic_link_exists;

#include <portcullis/intron.hpp>
#include <portcullis/portcullis_fs.hpp>
#include <portcullis/junction_system.hpp>
#include <portcullis/performance.hpp>
#include <portcullis/rule_parser.hpp>
using portcullis::PortcullisFS;
using portcullis::Intron;
using portcullis::IntronHasher;
using portcullis::Performance;
using portcullis::JuncResultMap;


namespace portcullis {
    
typedef boost::error_info<struct JuncFilterError,string> JuncFilterErrorInfo;
struct JuncFilterException: virtual boost::exception, virtual std::exception { };


const string DEFAULT_FILTER_OUTPUT = "portcullis_filter/portcullis";
const string DEFAULT_FILTER_SOURCE = "portcullis";
const string DEFAULT_FILTER_RULE_FILE = "default_filter.json";
const string DEFAULT_FILTER_MODEL_FILE = "default_model.forest";
const string ST_IPOS_RULES_FILE = "selftrain_initial_pos.json";
const string ST_INEG_RULES_FILE = "selftrain_initial_neg.json";
const uint16_t DEFAULT_FILTER_THREADS = 1;


class JunctionFilter {

private:
    
    path junctionFile;
    path modelFile;
    path filterFile;
    path genuineFile;
    path referenceFile;
    path output;
    bool train;
    uint16_t threads;
    bool saveBad;
    int32_t maxLength;
    bool filterCanonical;
    bool filterSemi;
    bool filterNovel;    
    string source;
    bool verbose;    
    
    
public:
    
    static path scriptsDir;
    static path dataDir;
    
    JunctionFilter( const path& _junctionFile, 
                    const path& _output);
    
    virtual ~JunctionFilter() {
    }
    
       

public:
    
    path getFilterFile() const {
        return filterFile;
    }

    void setFilterFile(path filterFile) {
        this->filterFile = filterFile;
    }

    path getJunctionFile() const {
        return junctionFile;
    }

    void setJunctionFile(path junctionFile) {
        this->junctionFile = junctionFile;
    }

    path getOutput() const {
        return output;
    }

    void setOutput(path output) {
        this->output = output;
    }
    
    path getGenuineFile() const {
        return genuineFile;
    }

    void setGenuineFile(path genuineFile) {
        this->genuineFile = genuineFile;
    }

    path getModelFile() const {
        return modelFile;
    }

    void setModelFile(path modelFile) {
        this->modelFile = modelFile;
    }

    path getReferenceFile() const {
        return referenceFile;
    }

    void setReferenceFile(path referenceFile) {
        this->referenceFile = referenceFile;
    }

    bool isTrain() const {
        return train;
    }

    void setTrain(bool train) {
        this->train = train;
    }

    bool isSaveBad() const {
        return saveBad;
    }

    void setSaveBad(bool saveBad) {
        this->saveBad = saveBad;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }
    
    string getSource() const {
        return source;
    }

    void setSource(string source) {
        this->source = source;
    }
    
    uint16_t getThreads() const {
        return threads;
    }

    void setThreads(uint16_t threads) {
        this->threads = threads;
    }
    
    void setCanonical(const string& canonical) {
        vector<string> modes;
        boost::split( modes, canonical, boost::is_any_of(","), boost::token_compress_on );

        if (modes.size() > 3) {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Canonical filter mode contains too many modes.  Max is 2.")));
        }
        
        if (modes.empty()) {
            this->filterCanonical = false;
            this->filterSemi = false;
            this->filterNovel = false;            
        }
        else {
            this->filterCanonical = true;
            this->filterSemi = true;
            this->filterNovel = true;

            for(auto& m : modes) {
                string n = boost::to_upper_copy(m);
                if (n == "OFF") {
                    this->filterCanonical = false;
                    this->filterSemi = false;
                    this->filterNovel = false;
                }
                else if (n == "C") {
                    this->filterCanonical = false;                    
                }
                else if (n == "S") {
                    this->filterSemi = false;
                }
                else if (n == "N") {
                    this->filterNovel = false;    
                }
            }
        }
    }

    bool doCanonicalFiltering() const {        
        return this->filterCanonical || this->filterSemi || this->filterNovel;
    }

    int32_t isMaxLength() const {
        return maxLength;
    }

    void setMaxLength(int32_t maxLength) {
        this->maxLength = maxLength;
    }
    
    path getIntitalPosRulesFile() const {
        return path(dataDir.string() + "/" + ST_IPOS_RULES_FILE);
    }

    path getIntitalNegRulesFile() const {
        return path(dataDir.string() + "/" + ST_INEG_RULES_FILE);
    }

    void filter();
    
 
protected:
    
    void forestPredict(const JunctionList& all, JunctionList& pass, JunctionList& fail, const uint32_t L95);

    shared_ptr<Performance> calcPerformance(const JunctionList& pass, const JunctionList& fail) {
        return calcPerformance(pass, fail, false);
    }
    shared_ptr<Performance> calcPerformance(const JunctionList& pass, const JunctionList& fail, bool invert);
    
    void printFilteringResults(const JunctionList& in, const JunctionList& pass, const JunctionList& fail, const string& prefix);
    
    void doRuleBasedFiltering(const path& ruleFile, const JunctionList& all, JunctionList& pass, JunctionList& fail, const string& prefix, JuncResultMap& resultMap);
    
    uint32_t calcIntronThreshold(const JunctionList& pass);
    
public:
  
    static string helpMessage() {
        return string("\nPortcullis Filter Mode Help.\n\n") +
                      "Filters out junctions that are unlikely to be genuine or that have too little\n" +
                      "supporting evidence.  The user can control three stages of the filtering\n" +
                      "process.  First the user can perform filtering based on a pre-defined random\n" +
                      "forest model.  Second the user can specify a configuration file described a\n" +
                      "set of filtering rules to apply.  Third, the user can directly through the\n" +
                      "command line filter based on junction (intron) length, or the canonical label.\n\n" +
                      "Usage: portcullis filter [options] <junction-file>\n\n" +
                      "Options";
    }
    
    static int main(int argc, char *argv[]);
};

}