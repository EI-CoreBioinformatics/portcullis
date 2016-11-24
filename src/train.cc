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
using std::shared_ptr;
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
#include <ranger/ForestRegression.h>

#include <portcullis/ml/performance.hpp>
#include <portcullis/ml/k_fold.hpp>
using portcullis::ml::Performance;
using portcullis::ml::PerformanceList;
using portcullis::ml::KFold;

#include <portcullis/junction_system.hpp>
using portcullis::JunctionSystem;
using portcullis::JunctionList;

#include "train.hpp"
using portcullis::Train;

portcullis::Train::Train(const path& _junctionFile, const path& _refFile) {
	junctionFile = _junctionFile;
	refFile = _refFile;
	outputPrefix = "";
	folds = DEFAULT_TRAIN_FOLDS;
	trees = DEFAULT_TRAIN_TREES;
	threads = DEFAULT_TRAIN_THREADS;
	fraction = DEFAULT_TRAIN_FRACTION;
	regressionMode = false;
	verbose = false;
}




void portcullis::Train::testInstance(shared_ptr<Forest> f, const JunctionList& x) {
	// Convert testing set junctions into feature vector
	ModelFeatures mf;
	Data* testingData = mf.juncs2FeatureVectors(x);
	vector<string> catvars;
	f->setPredictionMode(true);
	f->setData(testingData, "Genuine", "", catvars);
	f->run(false);
	delete testingData;
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
	for (size_t i = 0; i < all_junctions.size(); i++) {
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
		cout << "Training on full dataset" << endl;
		ForestPtr f = ModelFeatures().trainInstance(junctions, JunctionList(), outputPrefix.string(), trees, threads, false, true, false, false);
		f->saveToFile();
		f->writeOutput(&cout);
	}
	// Assess performance of the model if requested
	// Makes no sense to do cross validation on less than 2-fold
	if (folds >= 2) {
		// Setup k fold cross validation to estimate real performance
		KFold<JunctionList::const_iterator> kf(folds, junctions.begin(), junctions.end());
		JunctionList test, train;
		PerformanceList perfs;
		cout << endl << "Starting " << folds << "-fold cross validation" << endl;
		std::ofstream resout(outputPrefix.string() + ".cv_results");
		cout << "Fold\t" << Performance::longHeader() << endl;
		resout << "Fold\t" << Performance::longHeader() << endl;
		for (uint16_t i = 1; i <= folds; i++) {
			cout << i << "\t";
			resout << i << "\t";
			cout.flush();
			resout.flush();
			// Populate train and test for this step
			kf.getFold(i, back_inserter(train), back_inserter(test));
			// Train on this particular set
			ForestPtr f = ModelFeatures().trainInstance(train, JunctionList(), outputPrefix.string(), trees, threads, false, false, false, false);
			// Test model instance
			testInstance(f, test);
			uint32_t tp = 0, tn = 0, fp = 0, fn = 0;
			for (size_t j = 0; j < test.size(); j++) {
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
			shared_ptr<Performance> p = make_shared<Performance>(tp, tn, fp, fn);
			cout << p->toLongString() << endl;
			resout << p->toLongString() << endl;
			perfs.add(p);
			// Clear the train and test vectors in preparation for the next step
			train.clear();
			test.clear();
		}
		cout << "Cross validation completed" << endl << endl;
		perfs.outputMeanPerformance(resout);
		resout.close();
		cout << endl << "Saved cross validation results to file " << outputPrefix.string() << ".cv_results" << endl;
	}
}

void portcullis::Train::getRandomSubset(const JunctionList& in, JunctionList& out) {
	// Create a list of randomly shuffled indices with maximum of "in".
	std::vector<unsigned int> indices(in.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::random_shuffle(indices.begin(), indices.end());
	// Calculate number required in output
	const uint32_t outSize = (uint32_t)((double)in.size() * fraction);
	// Populate output
	for (size_t i = 0; i < outSize; i++) {
		out.push_back(in[indices[i]]);
	}
}

int portcullis::Train::main(int argc, char *argv[]) {
	// Portcullis args
	path junctionFile;
	path output;
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
	po::options_description generic_options(helpMessage(), w.ws_col, (unsigned)((double)w.ws_col / 1.7));
	generic_options.add_options()
	("output,o", po::value<path>(&output)->default_value(DEFAULT_TRAIN_OUTPUT),
	 "File name prefix for the random forest produced by this tool.")
	("reference,r", po::value<path>(&refFile),
	 "Either a reference bed file containing genuine junctions or file containing a line separated list of 1/0 corresponding to each entry in the input junction file indicating whether that entry is or isn't a genuine junction")
	("folds,k", po::value<uint16_t>(&folds)->default_value(DEFAULT_TRAIN_FOLDS),
	 "The level of cross validation to perform.  A value of 1 or less means do not do cross validation.  The default level of 5 is sufficient to get a reasonable feel for the accuracy of the model on portcullis datasets.")
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
	trainer.setOutputPrefix(output);
	trainer.setFolds(folds);
	trainer.setTrees(trees);
	trainer.setThreads(threads);
	trainer.setFraction(fraction);
	trainer.setVerbose(verbose);
	trainer.train();
	return 0;
}




