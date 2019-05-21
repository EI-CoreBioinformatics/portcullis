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

#include <portcullis/bam/genome_mapper.hpp>
using portcullis::bam::GenomeMapper;

#include <portcullis/ml/performance.hpp>
#include <portcullis/ml/model_features.hpp>
using portcullis::ml::Performance;
using portcullis::ml::ModelFeatures;

#include <portcullis/intron.hpp>
#include <portcullis/portcullis_fs.hpp>
#include <portcullis/junction_system.hpp>
using portcullis::PortcullisFS;
using portcullis::Intron;
using portcullis::IntronHasher;

#include "prepare.hpp"
using portcullis::PreparedFiles;


namespace portcullis {

    typedef boost::error_info<struct JuncFilterError, string> JuncFilterErrorInfo;

    struct JuncFilterException : virtual boost::exception, virtual std::exception {
    };


    const string DEFAULT_FILTER_OUTPUT = "portcullis_filter/portcullis";
    const string DEFAULT_FILTER_SOURCE = "portcullis";
    const string DEFAULT_FILTER_RULE_FILE = "default_filter.json";
    const string DEFAULT_FILTER_MODEL_FILE = "default_model.forest";
    const string ST_IPOS_RULES_FILE = "selftrain_initial_pos";
    const string ST_INEG_RULES_FILE = "selftrain_initial_neg";
    const uint16_t DEFAULT_FILTER_THREADS = 1;
    const uint16_t DEFAULT_SELFTRAIN_TREES = 250;
    const double DEFAULT_FILTER_THRESHOLD = 0.5;

    class JunctionFilter {
    private:

        path junctionFile;
        PreparedFiles prepData;
        path modelFile;
        path filterFile;
        path genuineFile;
        path referenceFile;
        path output;
        bool train;
        uint16_t threads;
        bool saveBad;
        bool saveFeatures;
        bool saveLayers;
        bool outputExonGFF;
        bool outputIntronGFF;
        uint32_t maxLength;
        bool filterCanonical;
        uint32_t minCov;
        bool filterSemi;
        bool filterNovel;
        string source;
        double threshold;
        bool smote;
        bool enn;
        bool precise;
        bool verbose;
        path initial;


    public:

        static path dataDir;

        JunctionFilter(const path& _prepDir,
                const path& _junctionFile,
                const path& _output,
                const path& _initial);

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

        double getThreshold() const {
            return threshold;
        }

        void setThreshold(double threshold) {
            this->threshold = threshold;
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

        bool doSaveBad() const {
            return saveBad;
        }

        void setSaveBad(bool saveBad) {
            this->saveBad = saveBad;
        }

        bool doSaveLayers() const {
            return saveLayers;
        }

        void setSaveLayers(bool saveLayers) {
            this->saveLayers = saveLayers;
        }

        bool isOutputExonGFF() const {
            return outputExonGFF;
        }

        void setOutputExonGFF(bool outputExonGFF) {
            this->outputExonGFF = outputExonGFF;
        }

        bool isOutputIntronGFF() const {
            return outputIntronGFF;
        }

        void setOutputIntronGFF(bool outputIntronGFF) {
            this->outputIntronGFF = outputIntronGFF;
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

        bool isENN() const {
            return enn;
        }

        void setENN(bool enn) {
            this->enn = enn;
        }

        bool isSmote() const {
            return smote;
        }

        void setSmote(bool smote) {
            this->smote = smote;
        }

        bool doSaveFeatures() const {
            return this->saveFeatures;
        }

        void setSaveFeatures(bool saveFeatures) {
            this->saveFeatures = saveFeatures;
        }

        void setCanonical(const string& canonical) {
            vector<string> modes;
            boost::split(modes, canonical, boost::is_any_of(","), boost::token_compress_on);
            if (modes.size() > 3) {
                BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                        "Canonical filter mode contains too many modes.  Max is 2.")));
            }
            if (modes.empty()) {
                this->filterCanonical = false;
                this->filterSemi = false;
                this->filterNovel = false;
            } else {
                this->filterCanonical = true;
                this->filterSemi = true;
                this->filterNovel = true;
                for (auto & m : modes) {
                    string n = boost::to_upper_copy(m);
                    if (n == "OFF") {
                        this->filterCanonical = false;
                        this->filterSemi = false;
                        this->filterNovel = false;
                    } else if (n == "C") {
                        this->filterCanonical = false;
                    } else if (n == "S") {
                        this->filterSemi = false;
                    } else if (n == "N") {
                        this->filterNovel = false;
                    }
                }
            }
        }

        bool doCanonicalFiltering() const {
            return this->filterCanonical || this->filterSemi || this->filterNovel;
        }

        uint32_t isMaxLength() const {
            return maxLength;
        }

        void setMaxLength(uint32_t maxLength) {
            this->maxLength = maxLength;
        }

        uint32_t getMinCov() const {
            return minCov;
        }

        void setMinCov(uint32_t minCov) {
            this->minCov = minCov;
        }

        bool isPrecise() const {
            return precise;
        }

        void setPrecise(bool precise) {
            this->precise = precise;
        }

        path getIntitalPosRulesFile(uint16_t index) const {
            return path(dataDir.string() + "/" + ST_IPOS_RULES_FILE + ".layer" + std::to_string(index) + ".json");
        }

        path getIntitalNegRulesFile(uint16_t index) const {
            return path(dataDir.string() + "/" + ST_INEG_RULES_FILE + ".layer" + std::to_string(index) + ".json");
        }

        void filter();


    protected:

        void forestPredict(const JunctionList& all, JunctionList& pass, JunctionList& fail, ModelFeatures& mf);

        shared_ptr<Performance> calcPerformance(const JunctionList& pass, const JunctionList& fail) {
            return calcPerformance(pass, fail, false);
        }
        shared_ptr<Performance> calcPerformance(const JunctionList& pass, const JunctionList& fail, bool invert);

        void printFilteringResults(const JunctionList& in, const JunctionList& pass, const JunctionList& fail, const string& prefix);

        void doRuleBasedFiltering(const path& ruleFile, const JunctionList& all, JunctionList& pass, JunctionList& fail);

        void categorise(shared_ptr<Forest> f, const JunctionList& all, JunctionList& pass, JunctionList& fail, double t);

        void createPositiveSet(const JunctionList& all, JunctionList& pos, JunctionList& unlabelled, ModelFeatures& mf);

        void createNegativeSet(uint32_t L95, const JunctionList& all, JunctionList& neg, JunctionList& failJuncs);

        double calcGoodThreshold(shared_ptr<Forest> f);

        void undersample(JunctionList& jl, size_t size);

	std::tuple<vector<string>, vector<string>> find_jsons(path ruleset);

	static bool sort_jsons(string& json1, string& json2);

    public:

        static string title() {
            return string("Portcullis Filter Mode Help");
        }

        static string description() {
            return string("Filters out junctions that are unlikely to be genuine or that have too little\n") +
                    "supporting evidence.  The user can control three stages of the filtering\n" +
                    "process.  First the user can perform filtering based on a random forest model\n" +
                    "self-trained on the provided data, alternatively the user can provide a pre-\n" +
                    "trained model.  Second the user can specify a configuration file describing a\n" +
                    "set of filtering rules to apply.  Third, the user can directly through the\n" +
                    "command line filter based on junction (intron) length, or the canonical label.\n\n" +
                    "This stage requires the prep directory and the tab file generated from the\n" +
                    "stage as input.";
        }

        static string usage() {
            return string("portcullis filter [options] <prep_data_dir> <junction_tab_file>");
        }

        static int main(int argc, char *argv[]);
    };

}
