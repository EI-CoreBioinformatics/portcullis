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
#include <boost/program_options.hpp>
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
namespace po = boost::program_options;
using boost::property_tree::ptree;
namespace qi    = boost::spirit::qi;
namespace phx   = boost::phoenix;

#include "intron.hpp"
#include "junction_system.hpp"
#include "portcullis_fs.hpp"
using portcullis::PortcullisFS;
using portcullis::Intron;
using portcullis::IntronHasher;

#include "junction_filter.hpp"


portcullis::Operator portcullis::stringToOp(const string& str) {
    if (boost::iequals(str, "EQ")) {
        return EQ;
    }
    else if (boost::iequals(str, "GT")) {
        return GT;
    }
    else if (boost::iequals(str, "LT")) {
        return LT;
    }
    else if (boost::iequals(str, "GTE")) {
        return GTE;
    }
    else if (boost::iequals(str, "LTE")) {
        return LTE;
    }
    else if (boost::iequals(str, "IN")) {
        return IN;
    }
    else if (boost::iequals(str, "NOT_IN")) {
        return NOT_IN;
    }
    else {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                        "Unrecognised operation: ") + str));
    }
}

string portcullis::opToString(const Operator op) {
    switch (op) {
        case EQ:
            return "EQ";                
        case GT:
            return "GT";
        case LT:
            return "LT";
        case GTE:
            return "GTE";
        case LTE:
            return "LTE";
        case IN:
            return "IN";
        case NOT_IN:
            return "NOT_IN";
        default:
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Unrecognised operation")));        
    }    
}

bool portcullis::isNumericOp(Operator op) {
    return (op != IN && op != NOT_IN);
}


portcullis::eval::eval(const NumericFilterMap& _numericmap, const SetFilterMap& _stringmap, const JunctionPtr _junc, JuncResultMap* _juncMap) : 
        boost::static_visitor<bool>() {
    numericmap = _numericmap;
    stringmap = _stringmap;
    junc = _junc;
    juncMap = _juncMap;
}

bool portcullis::eval::operator()(const var& v) const { 

    if (v=="T" || v=="t" || v=="true" || v=="True")
        return true;
    else if (v=="F" || v=="f" || v=="false" || v=="False")
        return false;
    else {
        // If it starts with an M then assume we are looking at a metric
        if (numericmap.count(v) > 0) {
            Operator op = numericmap.at(v).first;
            double threshold = numericmap.at(v).second;
            double value = getNumericFromJunc(v);
            bool res = evalNumberLeaf(op, threshold, value);
            if (!res) {
                juncMap->at(*(junc->getIntron())).push_back(v + " " + opToString(op) + " " + lexical_cast<string>(threshold));
            }
            return res;
        }
        else if (stringmap.count(v) > 0) {
            Operator op = stringmap.at(v).first;
            unordered_set<string> set = stringmap.at(v).second;

            string setstring = boost::algorithm::join(set, ", ");

            string value = getStringFromJunc(v);
            bool res = evalSetLeaf(op, set, value);
            if (!res) {
                juncMap->at(*(junc->getIntron())).push_back(v + " " + opToString(op) + " " + setstring);
            }
            return res;
        }
        else {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                "Unrecognised param: ") + v));
        }

    }
    return boost::lexical_cast<bool>(v); 
}

    
double portcullis::eval::getNumericFromJunc(const var& fullname) const {

    size_t pos = fullname.find(".");

    string name = pos == string::npos ? fullname : fullname.substr(0, pos);

    uint16_t metric_index = -1;
    for(size_t i = 0; i < NB_METRICS; i++) {
        if (boost::iequals(name, METRIC_NAMES[i])) {
            metric_index = i;
        }
    }

    if (metric_index == -1) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                "Unrecognised metric: ") + name));
    }

    switch(metric_index) {
        case 0:
            return 0.0;
        case 1:
            return (double)junc->getNbJunctionAlignments();
        case 2:
            return (double)junc->getNbDistinctAlignments();
        case 3:
            return (double)junc->getNbReliableAlignments();
        case 4:
            return (double)junc->getIntronSize();
        case 5:
            return (double)junc->getLeftAnchorSize();
        case 6:
            return (double)junc->getRightAnchorSize();
        case 7:
            return (double)junc->getMaxMinAnchor();
        case 8:
            return (double)junc->getDiffAnchor();
        case 9:
            return (double)junc->getNbDistinctAnchors();
        case 10:
            return junc->getEntropy();
        case 11:
            return (double)junc->getMaxMMES();
        case 12:
            return (double)junc->getHammingDistance5p();
        case 13:
            return (double)junc->getHammingDistance3p();
        case 14:
            return junc->getCoverage();
        case 15:
            return junc->isUniqueJunction() ? 1.0 : 0.0;
        case 16:
            return junc->isPrimaryJunction() ? 1.0 : 0.0;
        case 17:
            return junc->getMultipleMappingScore();
        case 18:
            return junc->getNbMismatches();
        case 19:
            return junc->getNbMultipleSplicedReads();
        case 20:
            return junc->getNbUpstreamJunctions();
        case 21:
            return junc->getNbDownstreamJunctions();
        case 22:
            return junc->getNbUpstreamFlankingAlignments();
        case 23:
            return junc->getNbDownstreamFlankingAlignments();

    }


    BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                        "Unrecognised metric")));        
}
    
string portcullis::eval::getStringFromJunc(const var& fullname) const {

    size_t pos = fullname.find(".");

    string name = pos == string::npos ? fullname : fullname.substr(0, pos);

    if (boost::iequals(name, "refname")) {
        return junc->getIntron()->ref.name;
    }
    else if (boost::iequals(name, "M1-canonical_ss")) {
        return string() + cssToChar(junc->getSpliceSiteType());
    }  

    BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                        "Unrecognised param: ") + name));        
}
    
bool portcullis::eval::evalNumberLeaf(Operator op, double threshold, double value) const {
    switch (op) {
        case EQ:
            return value == threshold;                
        case GT:
            return value > threshold;
        case LT:
            return value < threshold;
        case GTE:
            return value >= threshold;
        case LTE:
            return value <= threshold;
        default:
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Unrecognised operation")));
    }
}
    
bool portcullis::eval::evalSetLeaf(Operator op, unordered_set<string>& set, string value) const {
    switch (op) {
        case IN:
            return set.find(value) != set.end();
        case NOT_IN:
            return set.find(value) == set.end();
        default:
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Unrecognised operation")));
    }
}



portcullis::JunctionFilter::JunctionFilter( const path& _junctionFile, 
                    const path& _filterFile, 
                    const path& _outputDir, 
                    const string& _outputPrefix, 
                    const bool _saveBad,
                    bool _verbose) {
    junctionFile = _junctionFile;
    filterFile = _filterFile;
    outputDir = _outputDir;
    outputPrefix = _outputPrefix;
    saveBad = _saveBad;
    verbose = _verbose;
}
    
    
void portcullis::JunctionFilter::filter() {

    // Test if provided genome exists
    if (!exists(junctionFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find junction file at: ") + junctionFile.string()));
    }

    // Test if provided filter config file exists
    if (!exists(filterFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find filter configuration file at: ") + filterFile.string()));
    }

    if (!exists(outputDir)) {
        if (!create_directories(outputDir)) {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not create output directory at: ") + outputDir.string()));
        }
    }        

    /*string tmpOutputTabFile = outputDir.string() + "/tmp.tab";

    string filterCmd = string("python3 ") + filterJuncsPy.string() +
           " --input " + junctionFile.string() + 
           " --json " + filterFile.string() +
           " --out " + tmpOutputTabFile;

    if (verbose) {
        filterCmd += " --debug";

        cout << "Applying filter script: " << filterCmd << endl;
    }


    int exitCode = system(filterCmd.c_str());                    

    if (exitCode != 0 || !exists(tmpOutputTabFile)) {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Failed to successfully filter: ") + junctionFile.string()));
    }*/

    cout << "Loading JSON config ...";

    ptree pt;
    boost::property_tree::read_json(filterFile.string(), pt);

    cout << " done." << endl << endl;

    cout << "Loading junctions ...";
    cout.flush();

    // Load junction system
    JunctionSystem originalJuncs(junctionFile);

    cout << " done." << endl << endl;

    // Record junctions that pass the filter
    JunctionSystem passJunc;
    JunctionSystem failJunc;

    cout << "Filtering junctions ...";
    cout.flush();

    NumericFilterMap numericFilters;
    SetFilterMap stringFilters;
    JuncResultMap junctionMap;

    for(ptree::value_type& v : pt.get_child("parameters")) {
        string name = v.first;
        Operator op = stringToOp(v.second.get_child("operator").data());
        if (isNumericOp(op)) {
            double threshold = lexical_cast<double>(v.second.get_child("value").data());     
            numericFilters[name] = pair<Operator,double>(op, threshold);
        }
        else {

            unordered_set<string> set;
            for (auto& item : v.second.get_child("value")) {
                string val = item.second.get_value<string>();
                set.insert(val); 
            }
            stringFilters[name] = pair<Operator,unordered_set<string>>(op, set);
        }
    }

    const string expression = pt.get_child("expression").data();

    for (JunctionPtr junc : originalJuncs.getJunctions()) {

        junctionMap[*(junc->getIntron())] = vector<string>();

        if (parse(expression, junc, numericFilters, stringFilters, &junctionMap)) {
            passJunc.addJunction(junc);
        }
        else {
            failJunc.addJunction(junc);
        }
    }        

    cout << " done." << endl;

    // Output stats
    size_t diff = originalJuncs.size() - passJunc.size();

    cout << endl 
         << "Original junction file contained " << originalJuncs.size() << " junctions." << endl
         << "Filtered out " << diff << " junctions." << endl 
         << passJunc.size() << " junctions remaining" << endl << endl
         << "Saving junctions passing filter to disk:" << endl;

    passJunc.saveAll(outputDir.string() + "/" + outputPrefix + ".pass");        

    if (saveBad) {
        cout << "Saving junctions failing filter to disk:" << endl;

        failJunc.saveAll(outputDir.string() + "/" + outputPrefix + ".fail");
    }

    saveResults(junctionMap);
}
    
void portcullis::JunctionFilter::saveResults(JuncResultMap& results) {

    // Print descriptive output to file
    ofstream out(outputDir.string() + "/" + outputPrefix + ".results");

    out << Intron::locationOutputHeader() << "\t" << "Filter results..." << endl;

    for(auto& kv: results) {

        Intron i = kv.first;

        out << i << "\t";

        if (kv.second.empty()) {
            out << "PASS";
        }
        else {
            out << boost::algorithm::join(kv.second, "\t");
        }
        out << endl;
    }

    out.close();
}
    
    
/**
 * This function evaluates the truth status of a row parameter given the configuration present in the JSON file.
 * @param op Operation to be considered
 * @param threshold Threshold
 * @param param Value
 * @return True if parameter passes operation and threshold, false otherwise
 */
bool portcullis::JunctionFilter::parse(const string& expression, JunctionPtr junc, NumericFilterMap& numericFilters, SetFilterMap& stringFilters, JuncResultMap* results) {

    typedef std::string::const_iterator it;
    it f(expression.begin()), l(expression.end());
    parser<it> p;

    expr result;
    bool ok = qi::phrase_parse(f,l,p,qi::space,result);

    if (!ok) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Invalid expression: ") + expression));
    }

    // Evaluate results
    return boost::apply_visitor(eval(numericFilters, stringFilters, junc, results), result);
}
    
        
int portcullis::JunctionFilter::main(int argc, char *argv[]) {

    // Portcullis args
    string junctionFile;
    string filterFile;

    string outputDir;
    string outputPrefix;
    bool saveBad;
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(helpMessage(), 100);
    generic_options.add_options()
            ("output_dir,o", po::value<string>(&outputDir)->default_value(DEFAULT_FILTER_OUTPUT_DIR), 
                "Output directory for files generated by this program.")
            ("output_prefix,p", po::value<string>(&outputPrefix)->default_value(DEFAULT_FILTER_OUTPUT_PREFIX), 
                "File name prefix for files generated by this program.")
            ("filter_file,f", po::value<string>(&filterFile)->default_value(JunctionFilter::defaultFilterFile.string()), 
                "The filter configuration file to use.")
            ("save_bad,b", po::bool_switch(&saveBad)->default_value(false),
                "Saves bad junctions (i.e. junctions that fail the filter), as well as good junctions (those that pass)")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("junction-file,g", po::value<string>(&junctionFile), "Path to the junction file to process.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("junction-file", 1);


    // Combine non-positional options
    po::options_description cmdline_options;
    cmdline_options.add(generic_options).add(hidden_options);

    // Parse command line
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    // Output help information the exit if requested
    if (help || argc <= 1) {
        cout << generic_options << endl;
        return 1;
    }



    auto_cpu_timer timer(1, "\nPortcullis junction filter completed.\nTotal runtime: %ws\n\n");        

    cout << "Running portcullis in junction filter mode" << endl
         << "--------------------------------" << endl << endl;

    // Create the prepare class
    JunctionFilter filter(junctionFile, filterFile, outputDir, outputPrefix, saveBad, verbose);
    filter.filter();

    return 0;
}

path portcullis::JunctionFilter::defaultFilterFile = DEFAULT_FILTER_FILE;
path portcullis::JunctionFilter::filterJuncsPy = "";
