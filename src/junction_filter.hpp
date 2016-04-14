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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/variant/recursive_wrapper.hpp>
#include <boost/lexical_cast.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using boost::filesystem::absolute;
using boost::filesystem::copy_file;
using boost::filesystem::remove;
using boost::filesystem::exists;
using boost::filesystem::create_symlink;
using boost::filesystem::create_directory;
using boost::filesystem::symbolic_link_exists;
using boost::property_tree::ptree;
namespace qi    = boost::spirit::qi;
namespace phx   = boost::phoenix;

#include <portcullis/intron.hpp>
#include <portcullis/portcullis_fs.hpp>
#include <portcullis/junction_system.hpp>
using portcullis::PortcullisFS;
using portcullis::Intron;
using portcullis::IntronHasher;


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


const unordered_set<string> numericParams = {
    "M2-nb_reads",
    "M3-nb_dist_aln",
    "M4-nb_rel_aln",
    "M5-intron_size",
    "M6-left_anc_size",
    "M7-right_anc_size",
    "M8-max_min_anc",
    "M9-dif_anc",
    "M10-dist_anc",
    "M11-entropy",
    "M12-maxmmes",
    "M13-hamming5p",
    "M14-hamming3p",
    "M15-coverage",
    "M16-uniq_junc",
    "M17-primary_junc",
    "M18-mm_score",
    "M19-nb_mismatches",
    "M20-nb_msrs",
    "M21-nb_up_juncs",
    "M22-nb_down_juncs",
    "M23-up_aln",
    "M24-down_aln"        
};

const unordered_set<string> stringParams = {
    "M1-canonical_ss",
    "refname"
};


struct op_or  {};
struct op_and {};
struct op_not {};

typedef std::string var; 
template <typename tag> struct binop;
template <typename tag> struct unop;

typedef boost::variant<var, 
        boost::recursive_wrapper<unop <op_not> >, 
        boost::recursive_wrapper<binop<op_and> >,
        boost::recursive_wrapper<binop<op_or> >
        > expr;

template <typename tag> struct binop
{
    explicit binop(const expr& l, const expr& r) : oper1(l), oper2(r) { }
    expr oper1, oper2;
};

template <typename tag> struct unop
{
    explicit unop(const expr& o) : oper1(o) { }
    expr oper1;
};

enum Operator {
    EQ,
    GT,
    LT,
    GTE,
    LTE,
    IN,
    NOT_IN
};

typedef unordered_map<string, pair<Operator, double>> NumericFilterMap;
typedef unordered_map<string, pair<Operator, unordered_set<string>>> SetFilterMap;
typedef map<Intron, vector<string>, IntronComparator> JuncResultMap;

Operator stringToOp(const string& str);

string opToString(const Operator op);

bool isNumericOp(Operator op);

struct eval : boost::static_visitor<bool> {
    
    eval(const NumericFilterMap& _numericmap, const SetFilterMap& _stringmap, const JunctionPtr _junc, JuncResultMap* _juncMap);

    //
    bool operator()(const var& v) const;

    bool operator()(const binop<op_and>& b) const {        
        bool op1Res = recurse(b.oper1);
        bool op2Res = recurse(b.oper2);
        return op1Res && op2Res;
    }
    bool operator()(const binop<op_or>& b) const {
        bool op1Res = recurse(b.oper1);
        bool op2Res = recurse(b.oper2);
        return op1Res || op2Res;
    }
    bool operator()(const unop<op_not>& u) const {
        return !recurse(u.oper1);
    } 
    
    double getNumericFromJunc(const var& fullname) const;
    
    string getStringFromJunc(const var& fullname) const;
    
    bool evalNumberLeaf(Operator op, double threshold, double value) const;
    
    bool evalSetLeaf(Operator op, unordered_set<string>& set, string value) const;

    private:
        
    NumericFilterMap numericmap;
    SetFilterMap stringmap;
    JunctionPtr junc;
    JuncResultMap* juncMap;
    
    template<typename T>
        bool recurse(T const& v) const 
        { return boost::apply_visitor(*this, v); }
};


template <typename It, typename Skipper = qi::space_type>
    struct parser : qi::grammar<It, expr(), Skipper>
{
        parser() : parser::base_type(expr_) {
            
            expr_  = or_.alias();

            or_  = (and_ >> '|'  >> or_ ) [ qi::_val = phx::construct<binop<op_or > >(qi::_1_type(), qi::_2_type()) ] | and_   [ qi::_val = qi::_1_type() ];
            and_ = (not_ >> '&' >> and_)  [ qi::_val = phx::construct<binop<op_and> >(qi::_1_type(), qi::_2_type()) ] | not_   [ qi::_val = qi::_1_type() ];
            not_ = ('!' > simple       )  [ qi::_val = phx::construct<unop <op_not> >(qi::_1_type())     ] | simple [ qi::_val = qi::_1_type() ];

            simple = (('(' > expr_ > ')') | var_);
            var_ = qi::lexeme[+(qi::alpha|qi::digit|qi::char_("-")|qi::char_("_")|qi::char_("."))];

            BOOST_SPIRIT_DEBUG_NODE(expr_);
            BOOST_SPIRIT_DEBUG_NODE(or_);
            BOOST_SPIRIT_DEBUG_NODE(and_);
            BOOST_SPIRIT_DEBUG_NODE(not_);
            BOOST_SPIRIT_DEBUG_NODE(simple);
            BOOST_SPIRIT_DEBUG_NODE(var_);
        }

        private:
        qi::rule<It, var() , Skipper> var_;
        qi::rule<It, expr(), Skipper> not_, and_, or_, simple, expr_; 
};




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
    
    void saveResults(const JunctionSystem& js, JuncResultMap& results);
    
 
protected:
    
    /**
     * This function evaluates the truth status of a row parameter given the configuration present in the JSON file.
     * @param op Operation to be considered
     * @param threshold Threshold
     * @param param Value
     * @return True if parameter passes operation and threshold, false otherwise
     */
    bool parse(const string& expression, JunctionPtr junc, NumericFilterMap& numericFilters, SetFilterMap& stringFilters, JuncResultMap* results);
        
    
    wchar_t* convertCharToWideChar(const char* c);
    
    void executePythonMLFilter(const path& mlOutputFile);
    
    void forestPredict(const JunctionList& all, JunctionList& pass, JunctionList& fail);
    
    void calcPerformance(const JunctionList& pass, const JunctionList& fail);
    
    void printFilteringResults(const JunctionList& in, const JunctionList& pass, const JunctionList& fail, const string& prefix);
    
    void doRuleBasedFiltering(const path& ruleFile, const JunctionList& all, JunctionList& pass, JunctionList& fail, const string& prefix, JuncResultMap& resultMap);
    
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