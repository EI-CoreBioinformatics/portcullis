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

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
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
using std::cout;
using std::cerr;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
namespace bfs = boost::filesystem;
namespace po = boost::program_options;


#include <ranger/ForestClassification.h>
#include <ranger/ForestRegression.h>
#include <ranger/DataDouble.h>

#include <portcullis/bam/genome_mapper.hpp>
#include <portcullis/intron.hpp>
#include <portcullis/junction.hpp>
#include <portcullis/junction_system.hpp>
#include <portcullis/portcullis_fs.hpp>
#include <portcullis/performance.hpp>
#include <portcullis/rule_parser.hpp>
using portcullis::PortcullisFS;
using portcullis::Intron;
using portcullis::IntronHasher;
using portcullis::Performance;
using portcullis::eval;
using portcullis::bam::GenomeMapper;

#include "junction_filter.hpp"

portcullis::JunctionFilter::JunctionFilter( const path& _junctionFile, 
                    const path& _output) {
    junctionFile = _junctionFile;
    genomeFile = "";
    modelFile = "";
    genuineFile = "";
    output = _output;
    filterFile = "";
    referenceFile = "";
    saveBad = false;
    threads = 1;
    maxLength = 0;
    filterCanonical = false;
    filterSemi = false;
    filterNovel = false;
    source = DEFAULT_FILTER_SOURCE;
    verbose = false;
    threshold = DEFAULT_FILTER_THRESHOLD;
}
    
    
void portcullis::JunctionFilter::filter() {

    path outputDir = output.parent_path();
    string outputPrefix = output.leaf().string();
    
    if (outputDir.empty()) {
        outputDir = ".";
    }
    
    // Test if provided genome exists
    if (!exists(junctionFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find junction file at: ") + junctionFile.string()));
    }

    // Test if provided filter config file exists
    if (!modelFile.empty() && !exists(modelFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find filter model file at: ") + modelFile.string()));
    }
    
    if (!filterFile.empty() && !exists(filterFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find filter configuration file at: ") + filterFile.string()));
    }
    
    if (!genuineFile.empty() && !exists(genuineFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find file containing marked junction labels at: ") + genuineFile.string()));
    }
    
    if (!referenceFile.empty() && !exists(referenceFile)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not find reference BED file at: ") + referenceFile.string()));
    }

    if (!exists(outputDir)) {
        if (!bfs::create_directories(outputDir)) {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Could not create output directory at: ") + outputDir.string()));
        }
    }
    else if (!bfs::is_directory(outputDir)) {
        BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "File exists with name of suggested output directory: ") + outputDir.string()));            
    }
    
    cout << "Loading junctions from " << junctionFile.string() << " ...";
    cout.flush();

    // Load junction system
    JunctionSystem originalJuncs(junctionFile);

    cout << " done." << endl
         << "Found " << originalJuncs.getJunctions().size() << " junctions." << endl << endl;

    unordered_set<string> ref;
    if (!referenceFile.empty()) {        
        ifstream ifs(referenceFile.c_str());

        string line;
        // Loop through until end of file or we move onto the next ref seq
        while ( std::getline(ifs, line) ) {
            boost::trim(line);

            vector<string> parts; // #2: Search for tokens
            boost::split( parts, line, boost::is_any_of("\t"), boost::token_compress_on );

            // Ignore any non-entry lines
            if (parts.size() == 12) {
                int end = std::stoi(parts[7]) - 1;  // -1 to get from BED to portcullis coords for end pos
                string key = parts[0] + "(" + parts[6] + "," + std::to_string(end) + ")" + parts[5];
                ref.insert(key);
            }
        }
    }
    
    vector<bool> genuine;
    if (!genuineFile.empty()) {
        
        cout << "Loading list of correct predictions of performance analysis ...";
        cout.flush();
        
        Performance::loadGenuine(genuineFile, genuine);
        
        if (genuine.size() != originalJuncs.getJunctions().size()) {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                    "Genuine file contains ") + lexical_cast<string>(genuine.size()) +
                    " entries.  Junction file contains " + lexical_cast<string>(originalJuncs.getJunctions().size()) + 
                    " junctions.  The number of entries in both files must be the same to assess performance."));
        }
        
        // Copy over results into junction list
        for(size_t i = 0; i < originalJuncs.getJunctions().size(); i++) {
            originalJuncs.getJunctionAt(i)->setGenuine(genuine[i]);
        }
        
        cout << " done." << endl << endl;
    }
    
    // Also keep the current list of junctions
    JunctionList currentJuncs;
    
    // Copy everything into passJunc to begin with
    for(auto& j : originalJuncs.getJunctions()) {
        currentJuncs.push_back(j);
    }
    
    // To be overridden if we are training
    ModelFeatures mf;
    mf.initGenomeMapper(this->genomeFile);
    
    
    if (train) {
        
        // The initial positive and negative sets
        JunctionList pos, neg;
         
        cout << "Self training mode activated." << endl << endl;
        
        createPositiveSet(currentJuncs, pos, mf);
        
        createNegativeSet(mf.L95, currentJuncs, neg);
        
        cout << "Initial training set consists of " << pos.size() << " positive and " << neg.size() << " negative junctions." << endl << endl;
        
        cout << "Training markov models ...";
        cout.flush();
        mf.trainCodingPotentialModel(pos);
        mf.trainSplicingModels(pos, neg);
        cout << " done." << endl << endl;
        
        
        
        // Build the training set by combining the positive and negative sets
        JunctionList training;
        training.reserve(pos.size() + neg.size());
        training.insert(training.end(), pos.begin(), pos.end());
        training.insert(training.end(), neg.begin(), neg.end());
        
        JunctionSystem trainingSystem(training);
        trainingSystem.sort();        
        
        cout << "Training Random Forest" << endl
             << "----------------------" << endl << endl;
        bool done = false;
        shared_ptr<Forest> forest = nullptr;
        while(!done) {
        
            forest = mf.trainInstance(trainingSystem.getJunctions(), output.string() + ".selftrain", DEFAULT_SELFTRAIN_TREES, threads, true, true);
            const vector<double> importance = forest->getVariableImportance();
            mf.resetActiveFeatureIndex();
            uint16_t min_importance_idx = 10000;
            double min_importance_val = 100000.0;
            uint16_t min_feature_idx = 10000;
            // This approach to feature selection (i.e. removing worst features by importance)
            // doesn't appear to have much impact on final results and is computationally
            // expensive.  Removing for now.
            /*for(size_t i = 0; i < importance.size(); i++) {
                double imp = importance[i];
                int16_t j = mf.getNextActiveFeatureIndex();
                if (imp < 0.1 && imp < min_importance_val) {
                    min_importance_idx = i;
                    min_importance_val = imp;
                    min_feature_idx = j;
                }
            }*/
            if (min_importance_idx != 10000) {
                mf.features[min_feature_idx].active = false;
                cout << "Deactivating feature: " << mf.features[min_feature_idx].name << " - " << min_importance_val << endl;                
            }
            else {
                done = true;
            }
        }
        
        const vector<double> importance = forest->getVariableImportance();
        bool foundIrrelevant = false;
        mf.resetActiveFeatureIndex();
        cout << "Active features remaining:" << endl;
        for(auto& i : importance) {
            int16_t j = mf.getNextActiveFeatureIndex();
            cout << mf.features[j].name << " - " << i << endl;
        }
        
        forest->saveToFile();
        forest->writeOutput(&cout);
        
        modelFile = output.string() + ".selftrain.forest";
        cout << endl;
    }
       
    // Manage a junction system of all discarded junctions
    JunctionSystem discardedJuncs;
    
    // Do ML based filtering if requested
    if(!modelFile.empty() && exists(modelFile)){
        cout << "Predicting valid junctions using random forest model" << endl 
             << "----------------------------------------------------" << endl << endl;
        
        JunctionList passJuncs;
        JunctionList failJuncs;
        
        forestPredict(currentJuncs, passJuncs, failJuncs, mf);
                
        printFilteringResults(currentJuncs, passJuncs, failJuncs, string("Random Forest filtering results"));
        
        // Reset currentJuncs
        currentJuncs.clear();
        for(auto& j : passJuncs) {
            currentJuncs.push_back(j);
        }
        
        for(auto& j : failJuncs) {
            discardedJuncs.addJunction(j);
        }        
    }
    
    
    // Do rule based filtering if requested
    if (!filterFile.empty() && exists(filterFile)) {        
        
        JunctionList passJuncs;
        JunctionList failJuncs;
        JuncResultMap resultMap;
        
        doRuleBasedFiltering(filterFile, currentJuncs, passJuncs, failJuncs, "Rule-based filtering", resultMap);
        
        RuleFilter::saveResults(path(output.string() + ".rule_filtering.results"), originalJuncs, resultMap);
        
        printFilteringResults(currentJuncs, passJuncs, failJuncs, string("Rule-based filtering"));        
        
        // Reset currentJuncs
        currentJuncs.clear();
        for(auto& j : passJuncs) {
            currentJuncs.push_back(j);
        }
        
        // Add to discarded
        for(auto& j : failJuncs) {
            discardedJuncs.addJunction(j);
        }        
    }
    
    if (maxLength > 0 || this->doCanonicalFiltering()) {
        
        JunctionList passJuncs;
        JunctionList failJuncs;
        
        for(auto& j : currentJuncs) {
            
            bool pass = true;
            if (maxLength > 0) {
                if (j->getIntronSize() > maxLength) {
                    pass = false;
                }
            }
            
            if (pass && this->doCanonicalFiltering()) {
                if (this->filterNovel && j->getSpliceSiteType() == CanonicalSS::NO) {
                    pass = false;
                }
                if (this->filterSemi && j->getSpliceSiteType() == CanonicalSS::SEMI_CANONICAL) {
                    pass = false;
                }
                if (this->filterCanonical && j->getSpliceSiteType() == CanonicalSS::CANONICAL) {
                    pass = false;
                }
            }

            if (pass) {
                passJuncs.push_back(j);
            }
            else {
                failJuncs.push_back(j);
                discardedJuncs.addJunction(j);
            }
        }
        
        printFilteringResults(currentJuncs, passJuncs, failJuncs, string("Post filtering (length and/or canonical) results"));
        
        // Reset currentJuncs
        currentJuncs.clear();
        for(auto& j : passJuncs) {
            currentJuncs.push_back(j);
        }
    }
    
    cout << endl;
    
    JunctionSystem filteredJuncs;
    JunctionSystem refKeptJuncs;
        
    if (currentJuncs.empty()) {
        cout << "WARNING: Filters discarded all junctions from input." << endl;
    }
    else {
    
        cout << "Recalculating junction grouping and distance stats based on new junction list that passed filters ...";
        cout.flush();

        for(auto& j : currentJuncs) {
            filteredJuncs.addJunction(j);
        }
        
        if (!referenceFile.empty()) {
            for(auto& j : discardedJuncs.getJunctions()) {
                if (ref.count(j->locationAsString()) > 0) {
                    filteredJuncs.addJunction(j);
                    refKeptJuncs.addJunction(j);
                }
            }
        }

        filteredJuncs.calcJunctionStats();

        cout << " done." << endl << endl;
        
        if (!referenceFile.empty()) {
            cout << "Brought back " << refKeptJuncs.size() << " junctions that were discarded by filters but were present in reference file." << endl << endl;
        }
    }
    
    printFilteringResults(  originalJuncs.getJunctions(), 
                            filteredJuncs.getJunctions(), 
                            discardedJuncs.getJunctions(), 
                            string("Overall results"));

    cout << endl << "Saving junctions passing filter to disk:" << endl;

    filteredJuncs.saveAll(outputDir.string() + "/" + outputPrefix + ".pass", source + "_pass");
    
    if (saveBad) {
        cout << "Saving junctions failing filter to disk:" << endl;

        discardedJuncs.saveAll(outputDir.string() + "/" + outputPrefix + ".fail", source + "_fail");
        
        if (!referenceFile.empty()) {
            cout << "Saving junctions failing filters but present in reference:" << endl;

            refKeptJuncs.saveAll(outputDir.string() + "/" + outputPrefix + ".ref", source + "_ref");
        }
    }
}

void portcullis::JunctionFilter::createPositiveSet(const JunctionList& all, JunctionList& pos, ModelFeatures& mf) {
    JuncResultMap resultMap;
         
    cout << "Creating initial positive set for training" << endl
         << "------------------------------------------" << endl << endl;

    if (!genuineFile.empty()) {
        cout << "Performance of each positive filter layer (Low FPs is important):" << endl;
        cout << "LAYER\t" << Performance::longHeader() << endl;
    } 
    JunctionList p1, p2, f1, f2;
    RuleFilter::filter(this->getIntitalPosRulesFile(1), all, p1, f1, "Creating initial positive set for training", resultMap);
    if (!genuineFile.empty()) {
        cout << "1\t" << calcPerformance(p1, f1, true)->toLongString() << endl;
    }
    RuleFilter::filter(this->getIntitalPosRulesFile(2), p1, p2, f2, "Creating initial positive set for training", resultMap);
    if (!genuineFile.empty()) {
        cout << "2\t" << calcPerformance(p2, f2, true)->toLongString() << endl;
    }

    const uint32_t L95 = mf.calcIntronThreshold(p2);
    const uint32_t L95x2 = L95 * 2;
    
    JunctionList passJuncs, failJuncs;
    for(auto& j : p2) {
        if (j->getIntronSize() <= L95x2) {
            passJuncs.push_back(j);
        }
        else {
            failJuncs.push_back(j);
        }
    }
    if (!genuineFile.empty()) {
        cout << "L95x2\t" << calcPerformance(passJuncs, failJuncs, true)->toLongString() << endl;
    }
        
    // Analyse positive set to get L0.05 of intron size
    cout << endl << "Length of intron at 95th percentile of positive set (L95): " << mf.calcIntronThreshold(passJuncs) << endl << endl;

    cout << "Saving initial positive set:" << endl;
    JunctionSystem isp(passJuncs);
    isp.saveAll(output.string() + ".selftrain.initialset.pos", "portcullis_isp");


    for(auto& j : passJuncs) {
        JunctionPtr copy = make_shared<Junction>(*j);
        copy->setGenuine(true);
        pos.push_back(copy);            
    }
    if (!genuineFile.empty()) {
        
        JunctionSystem invalidPos;
        for(auto& j : passJuncs) {
            if (!j->isGenuine()) {
                invalidPos.addJunction(j);
            }
        }
        JunctionSystem missedPos;
        for(auto& j : failJuncs) {
            if (j->isGenuine()) {
                missedPos.addJunction(j);
            }
        }

        cout << "Saving invalid junctions in initial positive set to disk:" << endl;
        invalidPos.saveAll(output.string() + ".selftrain.initialset.invalidpos", "portcullis_invalid_isp");

        cout << "Saving missed positive junctions to disk:" << endl;
        missedPos.saveAll(output.string() + ".selftrain.initialset.missedpos", "portcullis_missed_isp");
    }
}

void portcullis::JunctionFilter::createNegativeSet(uint32_t L95, const JunctionList& all, JunctionList& neg) {
    
    JuncResultMap resultMap;
        
    cout << "Creating initial negative set for training" << endl
         << "------------------------------------------" << endl << endl;

    if (!genuineFile.empty()) {
       cout << "Performance of each negative filter layer (Low FNs is important):" << endl;
       cout << "LAYER\t" << Performance::longHeader() << endl;
    }
    JunctionList p1, p2, p3, p4, p5, p6, p7, f1, f2, f3, f4, f5, f6, f7;
    RuleFilter::filter(this->getIntitalNegRulesFile(1), all, p1, f1, "Creating initial negative set for training", resultMap);
    if (!genuineFile.empty()) {
       cout << "1\t" << calcPerformance(p1, f1, true)->toLongString() << endl;
    }
    RuleFilter::filter(this->getIntitalNegRulesFile(2), f1, p2, f2, "Creating initial negative set for training", resultMap);
    if (!genuineFile.empty()) {
       cout << "2\t" << calcPerformance(p2, f2, true)->toLongString() << endl;
    }
    RuleFilter::filter(this->getIntitalNegRulesFile(3), f2, p3, f3, "Creating initial negative set for training", resultMap);
    if (!genuineFile.empty()) {
       cout << "3\t" << calcPerformance(p3, f3, true)->toLongString() << endl;
    }
    RuleFilter::filter(this->getIntitalNegRulesFile(4), f3, p4, f4, "Creating initial negative set for training", resultMap);
    if (!genuineFile.empty()) {
       cout << "4\t" << calcPerformance(p4, f4, true)->toLongString() << endl;
    }
    RuleFilter::filter(this->getIntitalNegRulesFile(5), f4, p5, f5, "Creating initial negative set for training", resultMap);
    if (!genuineFile.empty()) {
       cout << "5\t" << calcPerformance(p5, f5, true)->toLongString() << endl;
    }
    RuleFilter::filter(this->getIntitalNegRulesFile(6), f5, p6, f6, "Creating initial negative set for training", resultMap);
    if (!genuineFile.empty()) {
       cout << "6\t" << calcPerformance(p6, f6, true)->toLongString() << endl;
    }
    
    JunctionList passJuncs, failJuncs;
    const uint32_t L95x5 = L95 * 5;
    for(auto& j : f6) {
       if (j->getIntronSize() > L95x5 && j->getMaxMMES() < 10 ) {
           p7.push_back(j);
       }
       else {
           failJuncs.push_back(j);
       }
    }
    if (!genuineFile.empty()) {
       cout << "L95x5\t" << calcPerformance(p7, failJuncs, true)->toLongString() << endl;
    }
    
    cout << endl << "Concatenating TNs from each layer to create negative set" << endl << endl;
    
    passJuncs.insert(passJuncs.end(), p1.begin(), p1.end());
    passJuncs.insert(passJuncs.end(), p2.begin(), p2.end());
    passJuncs.insert(passJuncs.end(), p3.begin(), p3.end());
    passJuncs.insert(passJuncs.end(), p4.begin(), p4.end());
    passJuncs.insert(passJuncs.end(), p5.begin(), p5.end());
    passJuncs.insert(passJuncs.end(), p6.begin(), p6.end());
    passJuncs.insert(passJuncs.end(), p7.begin(), p7.end());

    if (!genuineFile.empty()) {
       cout << "Final\t" << calcPerformance(passJuncs, failJuncs, true)->toLongString() << endl;
    }
    
    cout << endl << "Saving initial negative set:" << endl;
    JunctionSystem isn(passJuncs);
    isn.saveAll(output.string() + ".selftrain.initialset.neg", "portcullis_isn");

    for(auto& j : passJuncs) {
       JunctionPtr copy = make_shared<Junction>(*j);
       copy->setGenuine(false);
       neg.push_back(copy);            
    }

    if (!genuineFile.empty()) {
    
       JunctionSystem invalidNeg;
       JunctionSystem missedNeg;

       for(auto& j : passJuncs) {
           if (j->isGenuine()) {
               invalidNeg.addJunction(j);
           }
       }
       for(auto& j : failJuncs) {
           if (!j->isGenuine()) {
               missedNeg.addJunction(j);
           }
       }

       cout << "Saving genuine valid junctions in initial negative set to disk:" << endl;
       invalidNeg.saveAll(output.string() + ".selftrain.initialset.invalidneg", "portcullis_invalid_isn");

       cout << "Saving missed negative junctions to disk:" << endl;
       missedNeg.saveAll(output.string() + ".selftrain.initialset.missedneg", "portcullis_missed_isn");
    }
}


void portcullis::JunctionFilter::printFilteringResults(const JunctionList& in, const JunctionList& pass, const JunctionList& fail, const string& prefix) {
    // Output stats
    size_t diff = in.size() - pass.size();

    cout << endl << prefix << endl
         << "-------------------------" << endl
         << "Input contained " << in.size() << " junctions." << endl
         << "Output contains " << pass.size() << " junctions." << endl
         << "Filtered out " << diff << " junctions." << endl;
    
    if (!genuineFile.empty() && exists(genuineFile)) {
        shared_ptr<Performance> p = calcPerformance(pass, fail);
        cout << Performance::longHeader() << endl;
        cout << p->toLongString() << endl << endl;
    }    
}

shared_ptr<Performance> portcullis::JunctionFilter::calcPerformance(const JunctionList& pass, const JunctionList& fail, bool invert) {
    
    uint32_t tp = 0, tn = 0, fp = 0, fn = 0;

    if (invert) {
        for(auto& j : pass) {
            if (!j->isGenuine()) tn++; else fn++;
        }

        for(auto& j : fail) {
            if (j->isGenuine()) tp++; else fp++;
        }        
    }
    else {
        for(auto& j : pass) {
            if (j->isGenuine()) tp++; else fp++;
        }

        for(auto& j : fail) {
            if (!j->isGenuine()) tn++; else fn++;
        }
    }
    
    return make_shared<Performance>(tp, tn, fp, fn);
}


void portcullis::JunctionFilter::doRuleBasedFiltering(const path& ruleFile, const JunctionList& all, JunctionList& pass, JunctionList& fail, const string& prefix, JuncResultMap& resultMap) {
    cout << "Loading JSON rule-based filtering config file: " << ruleFile.string() << endl;

    cout << "Filtering junctions ...";
    cout.flush();

    map<string,int> filterCounts = RuleFilter::filter(ruleFile, all, pass, fail, prefix, resultMap);

    cout << " done." << endl << endl
         << "Number of junctions failing for each filter: " << endl;

    for(map<string,int>::iterator iterator = filterCounts.begin(); iterator != filterCounts.end(); iterator++) {        
        cout << iterator->first << ": " << iterator->second << endl;
    }

}
    
void portcullis::JunctionFilter::forestPredict(const JunctionList& all, JunctionList& pass, JunctionList& fail, ModelFeatures& mf) {
    
    cout << "Creating feature vector" << endl;
    
    Data* testingData = mf.juncs2FeatureVectors(all);
    
    cout << "Initialising random forest" << endl;
    
    shared_ptr<Forest> f = nullptr;
    if (train) {
        f = make_shared<ForestRegression>();
    }
    else {
        f = make_shared<ForestClassification>();
    }
    
    vector<string> catVars;
    
    f->init(
        "Genuine",                  // Dependant variable name
        MEM_DOUBLE,                 // Memory mode
        testingData,                // Data object
        0,                          // M Try (0 == use default)
        "",                         // Output prefix 
        DEFAULT_SELFTRAIN_TREES,                        // Number of trees (will be overwritten when loading the model)
        1234567890,                 // Seed for random generator               
        threads,                    // Number of threads
        IMP_GINI,                   // Importance measure 
        train ? DEFAULT_MIN_NODE_SIZE_REGRESSION : DEFAULT_MIN_NODE_SIZE_CLASSIFICATION,  // Min node size
        "",                         // Status var name 
        true,                       // Prediction mode
        true,                       // Replace 
        catVars,                    // Unordered categorical variable names (vector<string>)
        false,                      // Memory saving
        DEFAULT_SPLITRULE,          // Split rule
        false,                      // predall
        1.0);                       // Sample fraction
    
    f->setVerboseOut(&cerr);
    
    // Load trees from saved model
    f->loadFromFile(modelFile.string());
    
    cout << "Making predictions" << endl;
    f->run(verbose);
    
    path scorepath = output.string() + ".scores";
    ofstream scoreStream(scorepath.c_str());
    
    scoreStream << "Score\t" << Intron::locationOutputHeader() << "\tStrand\tSS\t" << testingData->getHeader() << endl;
    
    for(size_t i = 0; i < all.size(); i++) {

        scoreStream << f->getPredictions()[i][0] << "\t" << *(all[i]->getIntron()) 
                    << "\t" << strandToChar(all[i]->getConsensusStrand())
                    << "\t" << cssToChar(all[i]->getSpliceSiteType())
                    << "\t" << testingData->getRow(i) << endl;
    }
    
    scoreStream.close();
    
    
    if (!genuineFile.empty() && exists(genuineFile)) {
        vector<double> thresholds;
        for(double i = 0.05; i <= 0.5; i += 0.01) {
            thresholds.push_back(i);
        }
        double max_mcc = 0.0;
        double max_f1 = 0.0;
        double best_t_mcc = 0.0;
        double best_t_f1 = 0.0;
        
        cout << "Threshold\t" << Performance::longHeader() << endl;
        for(auto& t : thresholds) {
            JunctionList pjl;
            JunctionList fjl;
            categorise(f, all, pjl, fjl, t);
            shared_ptr<Performance> perf = calcPerformance(pjl, fjl);
            double mcc = perf->getMCC();
            double f1 = perf->getF1Score();
            cout << t << "\t" << perf->toLongString() << endl;
            if (mcc > max_mcc) {
                max_mcc = mcc;
                best_t_mcc = t;
            }
            if (f1 > max_f1) {
                max_f1 = f1;
                best_t_f1 = t;
            }
        }
        
        cout << "The best F1 score of " << max_f1 << " is achieved with threshold set at " << best_t_f1 << endl;
        cout << "The best MCC score of " << max_mcc << " is achieved with threshold set at " << best_t_mcc << endl;
        
        categorise(f, all, pass, fail, best_t_mcc);
    }
    else {
        categorise(f, all, pass, fail, threshold);
    }
    
    cout << "Saved junction scores to: " << scorepath << endl;
    
    
    delete testingData;
}

void portcullis::JunctionFilter::categorise(shared_ptr<Forest> f, const JunctionList& all, JunctionList& pass, JunctionList& fail, double t) {
    
    for(size_t i = 0; i < all.size(); i++) {
        if (f->getPredictions()[i][0] >= t) {
            pass.push_back(all[i]);
        }
        else {
            fail.push_back(all[i]);
        }
    }
}

int portcullis::JunctionFilter::main(int argc, char *argv[]) {

    path junctionFile;
    path genomeFile;
    path modelFile;
    path genuineFile;
    path filterFile;
    path referenceFile;
    path output;
    uint16_t threads;
    bool no_ml;
    bool saveBad;
    int32_t max_length;
    string canonical;
    string source;
    double threshold;
    bool verbose;
    bool help;
    
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);


    // Declare the supported options.
    po::options_description generic_options(helpMessage(), w.ws_col, w.ws_col/1.7);
    generic_options.add_options()
            ("output,o", po::value<path>(&output)->default_value(DEFAULT_FILTER_OUTPUT), 
                "Output prefix for files generated by this program.")
            ("filter_file,f", po::value<path>(&filterFile), 
                "If you wish to custom rule-based filter the junctions file, use this option to provide a list of the rules you wish to use.  By default we don't filter using a rule-based method, we instead filter via a self-trained random forest model.  See manual for more details.")
            ("model_file,m", po::value<path>(&modelFile), 
                "If you wish to use a custom random forest model to filter the junctions file, rather than self-training on the input dataset use this option to. See manual for more details.")
            ("genuine,g", po::value<path>(&genuineFile),
                "If you have a list of line separated boolean values in a file, indicating whether each junction in your input is genuine or not, then we can use that information here to gauge the accuracy of the predictions.")
            ("reference,r", po::value<path>(&referenceFile),
                "Reference annotation of junctions in BED format.  Any junctions found by the junction analysis tool will be preserved if found in this reference file regardless of any other filtering criteria.  If you need to convert a reference annotation from GTF or GFF to BED format portcullis contains scripts for this.")
            ("no_ml,n", po::bool_switch(&no_ml)->default_value(false),
                "Disables machine learning filtering")
            ("save_bad,b", po::bool_switch(&saveBad)->default_value(false),
                "Saves bad junctions (i.e. junctions that fail the filter), as well as good junctions (those that pass)")
            ("source", po::value<string>(&source)->default_value(DEFAULT_FILTER_SOURCE),
                "The value to enter into the \"source\" field in GFF files.")
            ("max_length,l", po::value<int32_t>(&max_length)->default_value(0),
                "Filter junctions longer than this value.  Default (0) is to not filter based on length.")
            ("canonical,c", po::value<string>(&canonical)->default_value("OFF"),
                "Keep junctions based on their splice site status.  Valid options: OFF,C,S,N. Where C = Canonical junctions (GU-AG), S = Semi-canonical junctions (AT-AC, or  GT-AG), N = Non-canonical.  OFF means, keep all junctions (i.e. don't filter by canonical status).  User can separate options by a comma to keep two categories.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_FILTER_THREADS), 
                "The number of threads to use during testing (only applies if using forest model).")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
        ("genome-file", po::value<path>(&genomeFile), "Path to the genome file to process.")        
        ("junction-file", po::value<path>(&junctionFile), "Path to the junction file to process.")
        ("threshold", po::value<double>(&threshold)->default_value(DEFAULT_FILTER_THRESHOLD), 
                "The threshold score at which we determine a junction to be genuine or not.")            
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("genome-file", 1);
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
         << "------------------------------------------" << endl << endl;

    // Create the prepare class
    JunctionFilter filter(junctionFile, output);
    filter.setSaveBad(saveBad);
    filter.setSource(source);
    filter.setVerbose(verbose);
    filter.setThreads(threads);
    filter.setMaxLength(max_length);
    filter.setCanonical(canonical);
    
    // Only set the filter rules if specified.
    filter.setFilterFile(filterFile);
    filter.setGenuineFile(genuineFile);
    if (modelFile.empty() && !no_ml) {
        filter.setTrain(true);    
    }
    else {
        filter.setTrain(false);
        if (!no_ml) {
            filter.setModelFile(modelFile);
        }
    }
    filter.setReferenceFile(referenceFile);
    filter.setGenomeFile(genomeFile);
    filter.setThreshold(threshold);
    
    filter.filter();

    return 0;
}

path portcullis::JunctionFilter::scriptsDir = ".";
path portcullis::JunctionFilter::dataDir = ".";
