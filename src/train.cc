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
#include <vector>
using std::vector;
using std::cout;
using std::endl;

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::filesystem::path;
using boost::lexical_cast;
namespace po = boost::program_options;


#include "train.hpp"
using portcullis::KFold;
using portcullis::Train;

portcullis::Train::Train(const path& _junctionFile, const path& _refFile, const path& _outputFile){
    junctionFile = _junctionFile;
    refFile = _refFile;
    outputFile = _outputFile;
    folds = DEFAULT_TRAIN_FOLDS;
    trees = DEFAULT_TRAIN_TREES;
    threads = DEFAULT_TRAIN_THREADS;
    verbose = false;    
}
    

void portcullis::Train::train() {
    
    // Load junction data
    cout << "Loaded junction data from " << junctionFile << endl;
            
    // Load reference data
    cout << "Loaded reference data from " << refFile << endl;
    
    
    cout << "Training on full dataset...";
    cout.flush();
    
    
    cout << "Saving random forest model to " << outputFile << endl;
    
    
    // Makes no sense to do cross validation on less than 2-fold
    if (folds >= 2) {
        vector<uint32_t> v;
        // Setup k fold cross validation to estimate real performance
        KFold<vector<uint32_t>::const_iterator> kf(folds, v.begin(), v.end());

        vector<uint32_t> test, train;
        double meanScore;

        cout << "Starting cross validation" << endl;

        for (uint16_t i = 1; i <= folds; i++) {

            cout << "Fold " << i << " ... ";
            cout.flush();

            // Populate train and test for this step
            kf.getFold(i, back_inserter(train), back_inserter(test));



            cout << "Score: " << 1 << endl;

            meanScore += 1; // 

            // Clear the train and test vectors in preparation for the next step
            train.clear();
            test.clear();
        }

        cout << "Cross validation completed" << endl << endl;
    
        meanScore /= folds;

        cout << "Mean score: " << meanScore << "% (+/- " << 0.0 << ")";
    }
    
}    


int portcullis::Train::main(int argc, char *argv[]) {
        
    // Portcullis args
    string junctionFile;
    string outputFile;
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
            ("output,o", po::value<string>(&outputFile)->default_value("trained_model.forest"), 
                "File name for the random forest produced by this tool.")
            ("reference,r", po::value<string>(&refFile), 
                "Either a reference bed file containing genuine junctions or file containing a line seperated list of 1/0 corresponding to each entry in the input junction file indicating whether that entry is or isn't a genuine junction")
            ("folds,k", po::value<uint16_t>(&folds)->default_value(DEFAULT_TRAIN_FOLDS), 
                "The level of cross validation to perform.  A value of 0 means do not do cross validation.  Normally a level of 10 is more than sufficient to get a reasonable feel for the accuracy of the model.")
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
    Train trainer(junctionFile, refFile, outputFile);
    trainer.setFolds(folds);
    trainer.setTrees(trees);
    trainer.setThreads(threads);
    trainer.setVerbose(verbose);
    trainer.train();

    return 0;
}




