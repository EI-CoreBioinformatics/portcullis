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
#include <iostream>
#include <memory>
#include <vector>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::make_shared;

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::filesystem::path;
using boost::lexical_cast;
namespace po = boost::program_options;

#include <ranger/globals.h>
#include <ranger/DataDouble.h>
#include <ranger/Forest.h>
#include <ranger/ForestClassification.h>

#include <portcullis/junction_system.hpp>
using portcullis::JunctionSystem;
using portcullis::JunctionList;

#include "train.hpp"
using portcullis::KFold;
using portcullis::Train;

portcullis::Train::Train(const path& _junctionFile, const path& _refFile){
    junctionFile = _junctionFile;
    refFile = _refFile;
    outputFile = "";
    folds = DEFAULT_TRAIN_FOLDS;
    trees = DEFAULT_TRAIN_TREES;
    threads = DEFAULT_TRAIN_THREADS;
    verbose = false;    
}
    
shared_ptr<Forest> portcullis::Train::trainInstance(const JunctionList& x) {
    
    // Convert junction list info to double*
    double* d = new double[x.size() * variableNames.size()];
    
    uint32_t row = 0;
    for (const auto& j : x) {        
        d[0 * x.size() + row] = (double)j->getNbJunctionAlignments();
        d[1 * x.size() + row] = (double)j->getNbDistinctAlignments();
        d[2 * x.size() + row] = (double)j->getNbReliableAlignments();
        d[3 * x.size() + row] = (double)j->getMaxMinAnchor();
        d[4 * x.size() + row] = (double)j->getDiffAnchor();
        d[5 * x.size() + row] = (double)j->getNbDistinctAnchors();
        d[6 * x.size() + row] = (double)j->getEntropy();
        d[7 * x.size() + row] = (double)j->getMaxMMES();
        d[8 * x.size() + row] = (double)j->getHammingDistance5p();
        d[9 * x.size() + row] = (double)j->getHammingDistance3p();
        d[10 * x.size() + row] = (double)j->isGenuine();
        
        row++;
    }
    
    Data* trainingData = new DataDouble(d, variableNames, x.size(), variableNames.size());
    
    shared_ptr<Forest> f = make_shared<ForestClassification>();
    
    vector<string> catVars;
    
    f->init(
        "Genuine",                  // Dependant variable name
        MEM_DOUBLE,                 // Memory mode
        trainingData,               // Data object
        0,                          // M Try
        outputFile.string(),        // Output prefix 
        trees,                      // Number of trees
        0,                          // Seed
        threads,                    // Number of threads
        DEFAULT_IMPORTANCE_MODE,    // Importance measure 
        0,                          // Target partition size
        "",                         // Status var name 
        false,                      // Prediction mode
        true,                       // Replace 
        catVars,                    // Unordered categorical variable names (vector<string>)
        false,                      // Memory saving
        DEFAULT_SPLITRULE,          // Split rule
        false,                      // predall
        1.0);                       // Sample fraction
            
    f->run(false);
    
    delete d;
    
    return f;
}

void portcullis::Train::testInstance(shared_ptr<Forest> f, const JunctionList& y) {
    
    // Convert junction list info to double*
    double* d = new double[y.size() * variableNames.size()];
    
    uint32_t row = 0;
    for (const auto& j : y) {        
        d[0 * y.size() + row] = (double)j->getNbJunctionAlignments();
        d[1 * y.size() + row] = (double)j->getNbDistinctAlignments();
        d[2 * y.size() + row] = (double)j->getNbReliableAlignments();
        d[3 * y.size() + row] = (double)j->getMaxMinAnchor();
        d[4 * y.size() + row] = (double)j->getDiffAnchor();
        d[5 * y.size() + row] = (double)j->getNbDistinctAnchors();
        d[6 * y.size() + row] = (double)j->getEntropy();
        d[7 * y.size() + row] = (double)j->getMaxMMES();
        d[8 * y.size() + row] = (double)j->getHammingDistance5p();
        d[9 * y.size() + row] = (double)j->getHammingDistance3p();
        d[10 * y.size() + row] = (double)j->isGenuine();
        
        row++;
    }
    
    Data* testingData = new DataDouble(d, variableNames, y.size(), variableNames.size());
    
    f->setPrediction_mode(true);   
    f->setData(testingData);
    f->run(false);
    
    delete d;
    
}

void portcullis::Train::train() {
    
    if (outputFile.empty() && folds < 2) {
        BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "You must specify either an output file to save the model too, and/or request a number of folds >= 2 for model evaluation")));
    }
    
    // Load junction data
    JunctionSystem junctions(junctionFile);
    cout << "Loaded " << junctions.size() << " junctions from " << junctionFile << endl;
           
    // Load reference data
    for(uint32_t i = 0; i < junctions.size(); i++) {
        junctions.getJunctionAt(i)->setGenuine(true);
    }
    cout << "Loaded reference data from " << refFile << endl;
    
    if (!outputFile.empty()) {
     
        cout << "Training on full dataset...";
        cout.flush();


        cout << "Saving random forest model to " << outputFile << endl;
    }
    
    // Assess performance of the model if requested
    // Makes no sense to do cross validation on less than 2-fold
    if (folds >= 2) {
        
        JunctionList jl = junctions.getJunctions();
        
        // Setup k fold cross validation to estimate real performance
        KFold<JunctionList::const_iterator> kf(folds, jl.begin(), jl.end());

        JunctionList test, train;
        vector<double> scores;

        cout << "Starting " << folds << "-fold cross validation" << endl;

        for (uint16_t i = 1; i <= folds; i++) {

            cout << "Fold " << i << " ... ";
            cout.flush();

            // Populate train and test for this step
            kf.getFold(i, back_inserter(train), back_inserter(test));
            
            cout << "Training size: " << train.size() << endl;
            cout << "Testing size: " << test.size() << endl;

            // Train on this particular set
            shared_ptr<Forest> f = trainInstance(train);
            
            f->writeOutput();
            
            // Test model instance
            for(const auto& j : test) {
                testInstance(f, test);
            }
            
            cout << "Score: " << 1 << endl;

            scores.push_back(1); // 

            // Clear the train and test vectors in preparation for the next step
            train.clear();
            test.clear();
        }

        cout << "Cross validation completed" << endl << endl;
    
        double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
        double mean = sum / scores.size();
        double sq_sum = std::inner_product(scores.begin(), scores.end(), scores.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / scores.size() - mean * mean);

        cout << "Mean score: " << mean << "% (+/- " << stdev << ")";
    }
    
}    


int portcullis::Train::main(int argc, char *argv[]) {
        
    // Portcullis args
    string junctionFile;
    string outputFile = "";
    string refFile;
    uint16_t folds;
    uint16_t trees;
    uint16_t threads;
    bool verbose;
    bool help;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    // Declare the supported options.
    po::options_description generic_options(helpMessage(), w.ws_col, (unsigned)((double)w.ws_col/1.7));
    generic_options.add_options()
            ("output,o", po::value<string>(&outputFile), 
                "File name for the random forest produced by this tool.")
            ("reference,r", po::value<string>(&refFile), 
                "Either a reference bed file containing genuine junctions or file containing a line separated list of 1/0 corresponding to each entry in the input junction file indicating whether that entry is or isn't a genuine junction")
            ("folds,k", po::value<uint16_t>(&folds)->default_value(DEFAULT_TRAIN_FOLDS), 
                "The level of cross validation to perform.  A value of 0 means do not do cross validation.  Normally a level of 5 is sufficient to get a reasonable feel for the accuracy of the model on portcullis datasets.")
            ("trees,n", po::value<uint16_t>(&trees)->default_value(DEFAULT_TRAIN_TREES), 
                "The number of trees to build in the random forest.  More trees will produce better results but at computational expense.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_TRAIN_THREADS), 
                "The number of threads to use during training.")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("junction-file", po::value<string>(&junctionFile), "Path to the junction file to process.")
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
    if (outputFile.empty()) {
        trainer.setOutputFile(outputFile);
    }    
    trainer.setFolds(folds);
    trainer.setTrees(trees);
    trainer.setThreads(threads);
    trainer.setVerbose(verbose);
    trainer.train();

    return 0;
}




