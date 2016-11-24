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
using std::cerr;
using std::endl;
using std::ostream;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>

#include <ranger/ForestProbability.h>

#include <portcullis/seq_utils.hpp>
using portcullis::SeqUtils;

#include <portcullis/ml/ss_forest.hpp>

portcullis::ml::SemiSupervisedForest::SemiSupervisedForest(ModelFeatures& mf,
		const JunctionList& _labelled, const JunctionList& _unlabelled,
		string _outputPrefix, uint16_t _trees, uint16_t _threads, double _contribution, bool _verbose) {
	verbose = _verbose;
	contribution = _contribution;
	labelled = mf.juncs2FeatureVectors(_labelled);
	if (verbose) cout << "Created labelled FV with " << labelled->getNumRows() << " entries." << endl;
	unlabelled = mf.juncs2FeatureVectors(_unlabelled);
	if (verbose) cout << "Created unlabelled FV with " << unlabelled->getNumRows() << " entries." << endl;
	// Combine labelled and unlabelled
	// This is inefficient but it makes life easier... hopefully it doesn't take too
	// much memory
	all = mf.juncs2FeatureVectors(_labelled, _unlabelled);
	if (verbose) cout << "Created combined FV with " << all->getNumRows() << " entries." << endl;
	outputPrefix = _outputPrefix;
	threads = _threads;
	trees = _trees;
}

portcullis::ml::SemiSupervisedForest::~SemiSupervisedForest() {
	delete labelled;
	delete unlabelled;
	delete all;
}

ForestPtr portcullis::ml::SemiSupervisedForest::train() {
	if (verbose) cout << "Initialising random forest" << endl;
	ForestPtr l = make_shared<ForestProbability>();
	vector<string> catVars;
	l->init(
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
		"",                         // Status var name (not required unless using survival)
		false,                      // Prediction mode
		true,                      // Sample with replacement
		catVars,                    // Unordered categorical variable names (vector<string>)
		false,                      // Memory saving
		DEFAULT_SPLITRULE,          // Split rule
		true,                       // Retain results for all trees (predict all)
		1.0);                       // Sample fraction
	vector<double> case_weights = {0.7, 0.3};
	//l->setCaseWeights(case_weights);
	l->setVerboseOut(&cout);
	if (verbose) cout << "Training on labelled data" << endl;
	l->run(verbose);
	string orig_forest = outputPrefix + ".ssrf.0.forest";
	l->saveToFile(orig_forest);
	vector<double> oobe;
	// Store out of box prediction error for first run on just labelled data
	oobe.push_back(l->getOverallPredictionError());
	cout << "OOBE: " << l->getOverallPredictionError() << endl;
	l->setPredictionMode(true);
	l->setData(unlabelled, "Genuine", "", catVars);
	l->run(verbose);
	ForestPtr u = make_shared<ForestProbability>();
	u->init(
		MEM_DOUBLE,                 // Memory mode
		0,                          // M Try (0 == use default)
		outputPrefix,               // Output prefix
		trees,                      // Number of trees
		1236456789,                 // Use fixed seed to avoid non-deterministic behaviour as much as possible
		threads,                    // Number of threads
		IMP_GINI,                   // Importance measure
		DEFAULT_MIN_NODE_SIZE_PROBABILITY,  // Min node size
		false,                      // Prediction mode
		true,                       // Replace
		false,                      // Memory saving
		DEFAULT_SPLITRULE,          // Split rule
		true,                       // predall
		1.0);                       // Sample fraction
	u->setData(all, "Genuine", "", catVars);
	// Loop until no improvement using deterministic annealing
	bool initOnly = false;
	bool first = true;
	bool improved = true;
	uint16_t repeat = 1;
	uint16_t it = 1;
	string best_forest;
	while (improved || repeat <= REPEAT_LIMIT) {
		// If OOBE improved during the last iteration make new predictions on
		// unlabelled data with the current model
		if (improved && !first) {
			if (verbose) cout << "Making predictions on the unlabelled set using current model" << endl;
			u->setPredictionMode(true);
			u->run(verbose);
			if (verbose) cout << "Made predictions." << endl;
		}
		// For the unlabelled set draw random labels using the actual probability
		// distributions of each tree in the forest and reassign genuine flag
		for (size_t i = 0; i < unlabelled->getNumRows(); i++) {
			// Note we are assuming field 0 contains the label
			// Override field 0 with the new label
			bool error = false;
			double pred = first ? l->makePrediction(i) : u->makePrediction(labelled->getNumRows() + i);
			//cout << pred << endl;
			all->set(0, labelled->getNumRows() + i, pred, error);
			if (error) {
				BOOST_THROW_EXCEPTION(SSRFException() << SSRFErrorInfo(string(
										  "Error setting label for initially unlabelled junction ") + std::to_string(i)));
			}
		}
		if (verbose) cout << "Re-training newly labelled data" << endl;
		u->setPredictionMode(false);
		u->run(verbose);
		cout << "OOBE: " << u->getOverallPredictionError() << endl;
		double error_delta = oobe.back() - u->getOverallPredictionError();
		if (error_delta <= 0.0 && !first) {
			improved = false;
			repeat++;
			cout << "No improvement with this permutation" << endl;
		}
		else {
			oobe.push_back(u->getOverallPredictionError());
			improved = true;
			repeat = 1;
			forest = u;
			if (!first) cout << "Improvement of " << error_delta << " to OOBE with this iteration" << endl;
			best_forest = outputPrefix + ".ssrf." + std::to_string(it++) + ".forest";
			u->saveToFile(best_forest);
		}
		first = false;
	}
	if (oobe.back() - oobe.front() >= 0.0) {
		initOnly = true;
		cout << "Using unlabelled data does not appear to have improved the model.  Reverting to original model derived only from initally labelled data." << endl;
	}
	else {
		for (size_t j = labelled->getNumRows(); j < labelled->getNumRows() + unlabelled->getNumRows(); j++) {
			bool error = false;
			all->set(0, j, all->get(j, 0), error);
		}
	}
	// Revert to the best forest
	ForestPtr b = make_shared<ForestProbability>();
	b->init(
		"Genuine",                  // Dependant variable name
		MEM_DOUBLE,                 // Memory mode
		initOnly ? labelled : all,                   // Data object
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
	b->run(verbose);
	return b;
}
