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

#include <sys/ioctl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
using std::ifstream;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::make_shared;

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;
namespace po = boost::program_options;

#include <ranger/globals.h>
#include <ranger/DataDouble.h>
#include <ranger/Forest.h>
#include <ranger/ForestClassification.h>

#include <portcullis/junction_system.hpp>
#include <portcullis/performance.hpp>
using portcullis::JunctionSystem;
using portcullis::JunctionList;
using portcullis::Performance;

#include "train.hpp"
using portcullis::KFold;
using portcullis::Train;

portcullis::Train::Train(const path& _junctionFile, const path& _refFile){
    junctionFile = _junctionFile;
    refFile = _refFile;
    outputPrefix = "";
    folds = DEFAULT_TRAIN_FOLDS;
    trees = DEFAULT_TRAIN_TREES;
    threads = DEFAULT_TRAIN_THREADS;
    fraction = DEFAULT_TRAIN_FRACTION;
    verbose = false;    
}
    
shared_ptr<Forest> portcullis::Train::trainInstance(const JunctionList& x) {
    
    // Convert junction list info to double*
    double* d = new double[x.size() * variableNames.size()];
    
    uint32_t row = 0;
    for (const auto& j : x) {        
        
        d[0 * x.size() + row] = j->getNbJunctionAlignments();
        d[1 * x.size() + row] = j->getNbDistinctAlignments();
        d[2 * x.size() + row] = j->getNbReliableAlignments();
        d[3 * x.size() + row] = j->getMaxMinAnchor();
        d[4 * x.size() + row] = j->getDiffAnchor();
        d[5 * x.size() + row] = j->getNbDistinctAnchors();
        d[6 * x.size() + row] = j->getEntropy();
        d[7 * x.size() + row] = j->getMaxMMES();
        d[8 * x.size() + row] = j->getHammingDistance5p();
        d[9 * x.size() + row] = j->getHammingDistance3p();
        d[10 * x.size() + row] = j->isGenuine();
        
        row++;
    }
    
    Data* trainingData = new DataDouble(d, variableNames, x.size(), variableNames.size());
    
    shared_ptr<Forest> f = make_shared<ForestClassification>();
    
    vector<string> catVars;
    
    f->init(
        "Genuine",                  // Dependant variable name
        MEM_DOUBLE,                 // Memory mode
        trainingData,               // Data object
        0,                          // M Try (0 == use default)
        outputPrefix.string(),      // Output prefix 
        trees,                      // Number of trees
        0,                          // Seed (0 == no random seed)
        threads,                    // Number of threads
        IMP_GINI,                   // Importance measure 
        DEFAULT_MIN_NODE_SIZE_CLASSIFICATION,  // Min node size
        "",                         // Status var name 
        false,                      // Prediction mode
        true,                       // Replace 
        catVars,                    // Unordered categorical variable names (vector<string>)
        false,                      // Memory saving
        DEFAULT_SPLITRULE,          // Split rule
        false,                      // predall
        1.0);                       // Sample fraction
            
    
    f->setVerboseOut(&cerr);
    f->run(false);
    
    delete d;
    
    return f;
}

void portcullis::Train::testInstance(shared_ptr<Forest> f, const JunctionList& x) {
    
    // Convert junction list info to double*
    double* d = new double[x.size() * variableNames.size()];
    
    uint32_t row = 0;
    for (const auto& j : x) {        
        d[0 * x.size() + row] = j->getNbJunctionAlignments();
        d[1 * x.size() + row] = j->getNbDistinctAlignments();
        d[2 * x.size() + row] = j->getNbReliableAlignments();
        d[3 * x.size() + row] = j->getMaxMinAnchor();
        d[4 * x.size() + row] = j->getDiffAnchor();
        d[5 * x.size() + row] = j->getNbDistinctAnchors();
        d[6 * x.size() + row] = j->getEntropy();
        d[7 * x.size() + row] = j->getMaxMMES();
        d[8 * x.size() + row] = j->getHammingDistance5p();
        d[9 * x.size() + row] = j->getHammingDistance3p();
        d[10 * x.size() + row] = j->isGenuine();
        
        row++;
    }
    
    Data* testingData = new DataDouble(d, variableNames, x.size(), variableNames.size());
    
    f->setPredictionMode(true);
    f->setData(testingData);
    f->run(false);
        
}

void portcullis::Train::train() {
    
    // Ensure output directory exists
    if (!outputPrefix.parent_path().empty()) {
        if (!bfs::exists(outputPrefix.parent_path())) {
            if (!bfs::create_directories(outputPrefix.parent_path())) {
                BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "Could not create output directory at: ") + outputPrefix.parent_path().string()));
            }
        }
    }
    
    if (outputPrefix.empty() && folds < 2) {
        BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "You must specify either an output file to save the model too, and/or request a number of folds >= 2 for model evaluation")));
    }
    
    if (!bfs::exists(junctionFile) && !bfs::symbolic_link_exists(junctionFile)) {
        BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "Junctions input file does not exist")));
    }
    
    if (!bfs::exists(refFile) && !bfs::symbolic_link_exists(refFile)) {
        BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "Reference file does not exist")));
    }
    
    if (fraction <= 0.0 || fraction > 1.0) {
        BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "Fraction: ") + lexical_cast<string>(fraction) + ". Valid values for \"fraction\" are between 0.0 and 1.0"));
    }
    
    
    // Load junction data
    JunctionSystem js;
    js.load(junctionFile, true);
    JunctionList all_junctions = js.getJunctions();
    cout << "Loaded " << all_junctions.size() << " junctions from " << junctionFile << endl;
           
    // Load reference data    
    vector<bool> genuine;
    Performance::loadGenuine(refFile, genuine);
    
    if (genuine.size() != all_junctions.size()) {
        BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "Ref data does not contain the same number of entries as the junctions input file.")));
    }
    
    // Copy over results into junction list
    for(size_t i = 0; i < all_junctions.size(); i++) {
        all_junctions[i]->setGenuine(genuine[i]);
    }
    
    cout << "Loaded reference data from " << refFile << endl;
    
    
    JunctionList junctions;
    
    if (fraction < 1.0) {
        
        getRandomSubset(all_junctions, junctions);
        
        cout << "Randomly subsampling input data to " << junctions.size() << " junctions." << endl;
    }
    else {
        junctions = all_junctions;
    }
    
    if (!outputPrefix.empty()) {
     
        cout << "Training on full dataset ...";
        cout.flush();

        shared_ptr<Forest> f = trainInstance(junctions);
        
        cout << " done." << endl;
        
        f->saveToFile();
        f->writeOutput();        
    }
    
    // Assess performance of the model if requested
    // Makes no sense to do cross validation on less than 2-fold
    if (folds >= 2) {
        
        // Setup k fold cross validation to estimate real performance
        KFold<JunctionList::const_iterator> kf(folds, junctions.begin(), junctions.end());

        JunctionList test, train;
        vector<double> f1s;
        vector<double> recalls;
        vector<double> precisions;

        cout << endl << "Starting " << folds << "-fold cross validation" << endl;
        cout << std::fixed << std::setprecision(2);

        for (uint16_t i = 1; i <= folds; i++) {

            cout << "Fold " << i << " ... ";
            cout.flush();

            // Populate train and test for this step
            kf.getFold(i, back_inserter(train), back_inserter(test));
            
            if (verbose) {
                cout << endl;
                cout << "Training size: " << train.size() << endl;
                cout << "Testing size: " << test.size() << endl;
            }

            // Train on this particular set
            shared_ptr<Forest> f = trainInstance(train);
            
            // Test model instance
            testInstance(f, test);
            
            uint32_t tp = 0, tn = 0, fp = 0, fn = 0;

            for(size_t j = 0; j < test.size(); j++) {
                double pred = f->getPredictions()[j][0];
                bool p = std::isnan(pred) ? false : pred == 1.0;
                bool r = test[j]->isGenuine();
                
                //cout << pred << " " << p << " " << r << endl;                 
                
                if (r) {
                    if (p) tp++; else fn++;
                }
                else {
                    if (p) fp++; else tn++;
                }
            }
            
            Performance p(tp, tn, fp, fn);
            
            
            cout << p.toString() << endl;

            f1s.push_back(p.getF1Score());
            recalls.push_back(p.getRecall());
            precisions.push_back(p.getPrecision());

            // Clear the train and test vectors in preparation for the next step
            train.clear();
            test.clear();
        }

        cout << "Cross validation completed" << endl << endl;
    
        outputMeanScore(recalls, "recall");
        outputMeanScore(precisions, "precisions");
        outputMeanScore(f1s, "F1");
    }
    
}

void portcullis::Train::outputMeanScore(const vector<double>& scores, const string& score_type) {
    double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
    double mean = sum / scores.size();
    double sq_sum = std::inner_product(scores.begin(), scores.end(), scores.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / scores.size() - mean * mean);

    cout << "Mean " << score_type << ": " << mean << "% (+/- " << stdev << "%)" << endl;
}

void portcullis::Train::getRandomSubset(const JunctionList& in, JunctionList& out) {

    // Create a list of randomly shuffled indices with maximum of "in".
    std::vector<unsigned int> indices(in.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    
    // Calculate number required in output
    const uint32_t outSize = (uint32_t)((double)in.size() * fraction);

    // Populate output
    for(size_t i = 0; i < outSize; i++) {
        out.push_back(in[indices[i]]);
    }    
}

int portcullis::Train::main(int argc, char *argv[]) {
        
    // Portcullis args
    path junctionFile;
    string outputPrefix = "";
    path refFile;
    uint16_t folds;
    uint16_t trees;
    uint16_t threads;
    double fraction;
    bool verbose;
    bool help;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    // Declare the supported options.
    po::options_description generic_options(helpMessage(), w.ws_col, (unsigned)((double)w.ws_col/1.7));
    generic_options.add_options()
            ("output,o", po::value<string>(&outputPrefix), 
                "File name prefix for the random forest produced by this tool.")
            ("reference,r", po::value<path>(&refFile), 
                "Either a reference bed file containing genuine junctions or file containing a line separated list of 1/0 corresponding to each entry in the input junction file indicating whether that entry is or isn't a genuine junction")
            ("folds,k", po::value<uint16_t>(&folds)->default_value(DEFAULT_TRAIN_FOLDS), 
                "The level of cross validation to perform.  A value of 0 means do not do cross validation.  Normally a level of 5 is sufficient to get a reasonable feel for the accuracy of the model on portcullis datasets.")
            ("trees,n", po::value<uint16_t>(&trees)->default_value(DEFAULT_TRAIN_TREES), 
                "The number of trees to build in the random forest.  More trees will produce better results but at computational expense.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_TRAIN_THREADS), 
                "The number of threads to use during training.")
            ("fraction,f", po::value<double>(&fraction)->default_value(1.0), 
                "Fraction of the input data to use of training.  Note this is NOT the training / testing set ratio.  This is essentially a method for subsampling the data.")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("junction-file", po::value<path>(&junctionFile), "Path to the junction file to process.")
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



    auto_cpu_timer timer(1, "\nPortcullis training completed.\nTotal runtime: %ws\n\n");        

    cout << "Running portcullis in training mode" << endl
         << "-----------------------------------" << endl << endl;

    // Create the prepare class
    Train trainer(junctionFile, refFile);
    if (!outputPrefix.empty()) {
        trainer.setOutputPrefix(outputPrefix);
    }    
    trainer.setFolds(folds);
    trainer.setTrees(trees);
    trainer.setThreads(threads);
    trainer.setFraction(fraction);
    trainer.setVerbose(verbose);
    trainer.train();

    return 0;
}




