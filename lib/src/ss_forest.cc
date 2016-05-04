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

#include <iostream>
using std::cout;
using std::endl;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>

#include <portcullis/seq_utils.hpp>
using portcullis::SeqUtils;

#include <portcullis/ml/ss_forest.hpp>

ForestPtr portcullis::ml::SSForest::train() {
    
    if (verbose) cout << "Initialising random forest" << endl;
    ForestPtr f = make_shared<ForestProbability>();
    
    vector<string> catVars;
    
    f->init(
        "Genuine",                  // Dependant variable name
        MEM_DOUBLE,                 // Memory mode
        labelled,                   // Data object
        0,                          // M Try (0 == use default)
        outputPrefix,               // Output prefix 
        trees,                      // Number of trees
        1236456789,                 // Use fixed seed to avoid non-deterministic behaviour as much as possible
        threads,                    // Number of threads
        IMP_GINI,                   // Importance measure 
        DEFAULT_MIN_NODE_SIZE_PROBABILITY,  // Min node size
        "",                         // Status var name 
        false,                      // Prediction mode
        true,                       // Replace 
        catVars,                    // Unordered categorical variable names (vector<string>)
        false,                      // Memory saving
        DEFAULT_SPLITRULE,          // Split rule
        true,                       // predall
        1.0);                       // Sample fraction
    
    f->setVerboseOut(&cerr);
    
    if (verbose) cout << "Training on labelled data" << endl;
    f->run(verbose);
    
    string best_forest = outputPrefix+"/ssrf.0.forest";    
    f->saveToFile(best_forest);        
    
    vector<double> oobe;
    
    // Store out of box prediction error for first run on just labelled data
    oobe.push_back(f->getOverallPredictionError());
    
    // Loop until no improvement using deterministic annealing
    bool improved = true;
    uint16_t repeat = 1;
    uint16_t it = 1;
    while(improved || repeat <= REPEAT_LIMIT) {
        
        // If OOBE improved during the last iteration make new predictions on 
        // unlabelled data with the current model
        if (improved) {
            if (verbose) cout << "Making predictions on the unlabelled set using current model" << endl;
            f->setData(unlabelled);        
            f->setPredictionMode(true);
            f->run(verbose);
            
            // For the unlabelled set draw random labels using the actual probability 
            // distributions of each tree in the forest and reassign genuine flag
            for(size_t i = 0; i < unlabelled.size(); i++) {
                double score = f->getPrediction(i);
                double genuine = f->generateNewPrediction(i);
                
            }
        }
        
        
        
        if (verbose) cout << "Re-training using labelled and unlabelled data" << endl;
        f->setPredictionMode(false);
        f->run(verbose);
        
        cout << "OOBE: " << f->getOverallPredictionError() << endl;
        
        double error_delta = oobe.back() - f->getOverallPredictionError();
        if (error_delta <= 0.0) {
            improved = false;
            repeat++;
            cout << "No improvement with this permutation" << endl;
        }
        else {
            oobe.push_back(f->getOverallPredictionError());
            improved = true;
            repeat = 1;
            forest = fm;
            cout << "Improvement of " << error_delta << " to OOBE with this iteration" << endl;
            best_forest = outputPrefix+"/ssrf." + std::to_string(it++) + ".forest";
            f->saveToFile(best_forest);
        }
    }
    
    // Revert to the best forest
    f->loadFromFile(best_forest);
    
    return f;
}
