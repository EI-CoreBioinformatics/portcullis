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
#include <random>
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
using std::to_string;
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
#include <ranger/ForestProbability.h>
#include <ranger/DataDouble.h>

#include <portcullis/intron.hpp>
#include <portcullis/junction.hpp>
#include <portcullis/junction_system.hpp>
#include <portcullis/portcullis_fs.hpp>
#include <portcullis/python_helper.hpp>
using portcullis::PortcullisFS;
using portcullis::Intron;
using portcullis::IntronHasher;
using portcullis::PyHelper;

#include "junction_filter.hpp"
#include "prepare.hpp"

portcullis::JunctionFilter::JunctionFilter(const path& _prepDir, const path& _junctionFile,
		const path& _output) {
	junctionFile = _junctionFile;
	prepData.setPrepDir(_prepDir);
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
	smote = true;
	enn = true;
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
	if (!bfs::exists(prepData.getGenomeFilePath())) {
		BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
								  "Could not find prepared genome file at: ") + prepData.getGenomeFilePath().string()));
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

    // Also keep a shortcut to the list of junctions
    JunctionList currentJuncs = originalJuncs.getJunctions();

	unordered_set<string> ref;
	if (!referenceFile.empty()) {
		cout << "Loading junctions from reference: " << referenceFile.string() << " ...";
		cout.flush();
		ifstream ifs(referenceFile.c_str());
		string line;
		// Loop through until end of file or we move onto the next ref seq
		while (std::getline(ifs, line)) {
			boost::trim(line);
			vector<string> parts; // #2: Search for tokens
			boost::split(parts, line, boost::is_any_of("\t"), boost::token_compress_on);
			// Ignore any non-entry lines
			if (parts.size() == 12) {
				int end = std::stoi(parts[7]) - 1; // -1 to get from BED to portcullis coords for end pos
				string key = parts[0] + "(" + parts[6] + "," + std::to_string(end) + ")" + parts[5];
				ref.insert(key);
			}
		}
		cout << " done." << endl
			 << "Found " << ref.size() << " junctions in reference." << endl << endl;
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
		for (size_t i = 0; i < originalJuncs.getJunctions().size(); i++) {
			originalJuncs.getJunctionAt(i)->setGenuine(genuine[i]);
		}
		cout << " done." << endl << endl;
	}

	// To be overridden if we are training
	ModelFeatures mf;
	mf.initGenomeMapper(prepData.getGenomeFilePath());
	mf.features[1].active = false; // NB USRS          (BAD)
	mf.features[2].active = false; // NB DISTRS        (BAD)
	//mf.features[3].active=false;      // NB RELRS         (GOOD)
	mf.features[4].active = false; // ENTROPY          (BAD - JO LOGDEV ARE BETTER)
	//mf.features[5].active = false;    // REL2RAW          (GOOD)
	mf.features[6].active = false; // MAXMINANC        (BAD - MAXMMES IS BETTER)
	//mf.features[7].active=false;      // MAXMMES          (GOOD)
	//mf.features[8].active=false;      // MEAN MISMATCH    (GOOD)
	//mf.features[9].active=false;      // INTRON           (GOOD)
	//mf.features[10].active=false;     // MIN_HAMM         (GOOD)
	mf.features[11].active = false; // CODING POTENTIAL (BAD)
	//mf.features[12].active=false;     // POS WEIGHTS      (GOOD)
	//mf.features[13].active=false;     // SPLICE SIGNAL    (GOOD)
	/*mf.features[14].active=false;     // JO LOGDEV FEATURES BETTER THAN ENTROPY
	mf.features[15].active=false;
	mf.features[16].active=false;
	mf.features[17].active=false;
	mf.features[18].active=false;
	mf.features[19].active=false;
	mf.features[20].active=false;
	mf.features[21].active=false;
	mf.features[22].active=false;
	mf.features[23].active=false;
	mf.features[24].active=false;
	mf.features[25].active=false;
	mf.features[26].active=false;
	mf.features[27].active=false;
	mf.features[28].active=false;
	mf.features[29].active=false;
	 */
	double ratio = 0.0;
	if (train) {
		// The initial positive and negative sets
        JunctionList unlabelled, unlabelled2;
		cout << "Self training mode activated." << endl << endl;

        path rf_script = path("portcullis") / "rule_filter.py";
        vector<string> args;
        args.push_back(rf_script.string());

        args.push_back("--pos_json");
        args.push_back(dataDir.string() + "/selftrain_initial_pos.layer1.json");
        args.push_back(dataDir.string() + "/selftrain_initial_pos.layer2.json");
        args.push_back(dataDir.string() + "/selftrain_initial_pos.layer3.json");

        args.push_back("--neg_json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer1.json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer2.json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer3.json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer4.json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer5.json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer6.json");
        args.push_back(dataDir.string() + "/selftrain_initial_neg.layer7.json");

        args.push_back("--prefix=" + output.string() + ".selftrain.initialset");
        args.push_back(junctionFile.string());

        char* char_args[50];

        for(size_t i = 0; i < args.size(); i++) {
            char_args[i] = strdup(args[i].c_str());
        }

        PyHelper::getInstance().execute(rf_script.string(), (int)args.size(), char_args);

        // Load junction system
        JunctionSystem posSystem(path(output.string() + ".selftrain.initialset.pos.junctions.tab"));
        JunctionSystem negSystem(path(output.string() + ".selftrain.initialset.neg.junctions.tab"));
        posSystem.sort();
        negSystem.sort();
        JunctionList pos = posSystem.getJunctions();
        JunctionList neg = negSystem.getJunctions();

        // Load positive and negative set from disk

		cout << "Initial training set consists of " << pos.size() << " positive and " << neg.size() << " negative junctions." << endl << endl;
		ratio = 1.0 - ((double) pos.size() / (double) (pos.size() + neg.size()));
		cout << "Pos to neg ratio: " << ratio << endl << endl;

        // Ensure the L95 is set to what the python script generated.
        std::ifstream isL95(output.string() + ".selftrain.initialset.L95_intron_size.txt");
        bool foundL95 = false;
        for (int i = 0; i < 2; i++) {
            string line;
            std::getline(isL95, line);
            if (i == 1) {
                mf.L95 = lexical_cast<uint32_t>(line);
                foundL95 = true;
                break;
            }
        }

        // Double check we got that correctly
        if (!foundL95) {
            BOOST_THROW_EXCEPTION(JuncFilterException() << JuncFilterErrorInfo(string(
                                  "Problem loading L95 value from disk: " + output.string() + ".L95_intron_size.txt")));
        }
        cout << "Confirming intron length L95 is: " << mf.L95 << endl;

        cout << "Feature learning from training set ...";
		cout.flush();
        mf.trainCodingPotentialModel(pos);
		mf.trainSplicingModels(pos, neg);
		cout << " done." << endl << endl;

		cout << "Training Random Forest" << endl
			 << "----------------------" << endl << endl;
        shared_ptr<Forest> forest = mf.trainInstance(pos, neg, output.string() + ".selftrain", DEFAULT_SELFTRAIN_TREES, threads, true, true, smote, enn, saveFeatures);
		forest->saveToFile();
		modelFile = output.string() + ".selftrain.forest";
		cout << endl;
	}
	// Manage a junction system of all discarded junctions
	JunctionSystem discardedJuncs;
	// Do ML based filtering if requested
	if (!modelFile.empty() && exists(modelFile)) {
		cout << "Predicting valid junctions using random forest model" << endl
			 << "----------------------------------------------------" << endl << endl;
		JunctionList passJuncs;
		JunctionList failJuncs;
		forestPredict(currentJuncs, passJuncs, failJuncs, mf);
		printFilteringResults(currentJuncs, passJuncs, failJuncs, string("Random Forest filtering results"));
		// Reset currentJuncs
		currentJuncs.clear();
		for (auto & j : passJuncs) {
			currentJuncs.push_back(j);
		}
		for (auto & j : failJuncs) {
			discardedJuncs.addJunction(j);
		}
	}

    if (currentJuncs.empty()) {
        cout << "WARNING: Trained model discarded all junctions from input.  Will not apply any further filters." << endl;
    }
    else {

        // Do rule based filtering if requested
        if (!filterFile.empty() && exists(filterFile)) {

            JunctionSystem remainingJuncs;
            for (auto & j : currentJuncs) {
                remainingJuncs.addJunction(j);
            }
            remainingJuncs.saveAll(output.string() + ".rules_in", source + "_rules", false, false, false);

            path rf_script = path("portcullis") / "rule_filter.py";
            vector<string> args;
            args.push_back(rf_script.string());

            args.push_back("--json=" + filterFile.string());
            args.push_back("--prefix=" + output.string() + ".rules_out");
            args.push_back("--save_failed");
            args.push_back(output.string() + ".rules_in.junctions.tab");

            char* char_args[50];

            for(size_t i = 0; i < args.size(); i++) {
                char_args[i] = strdup(args[i].c_str());
            }

            PyHelper::getInstance().execute(rf_script.string(), (int)args.size(), char_args);

            // Load junction system
            JunctionSystem posSystem(path(output.string() + ".rules_out.passed.junctions.tab"));
            JunctionSystem negSystem(path(output.string() + ".rules_out.failed.junctions.tab"));
            posSystem.sort();
            negSystem.sort();

            // Reset currentJuncs
            currentJuncs.clear();
            for (auto & j : posSystem.getJunctions()) {
                currentJuncs.push_back(j);
            }
            // Add to discarded
            for (auto & j : negSystem.getJunctions()) {
                discardedJuncs.addJunction(j);
            }
        }

        if (currentJuncs.empty()) {
            cout << "WARNING: Rule-based filter discarded all junctions from input.  Will not apply any further filters." << endl;
        }
        else {

            if (maxLength > 0 || this->doCanonicalFiltering() || minCov > 1) {
                JunctionList passJuncs;
                JunctionList failJuncs;
                for (auto & j : currentJuncs) {
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
                    if (pass && this->getMinCov() > j->getNbSplicedAlignments()) {
                        pass = false;
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
                for (auto & j : passJuncs) {
                    currentJuncs.push_back(j);
                }
            }
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
		for (auto & j : currentJuncs) {
			filteredJuncs.addJunction(j);
		}
		uint32_t inref = 0;
		if (!referenceFile.empty()) {
			for (auto & j : currentJuncs) {
				if (ref.count(j->locationAsString()) > 0) {
					inref++;
				}
			}
			for (auto & j : discardedJuncs.getJunctions()) {
				if (ref.count(j->locationAsString()) > 0) {
					filteredJuncs.addJunction(j);
					refKeptJuncs.addJunction(j);
					inref++;
				}
			}
		}
		filteredJuncs.calcJunctionStats();
		cout << " done." << endl << endl;
		if (!referenceFile.empty()) {
			cout << "Brought back " << refKeptJuncs.size() << " junctions that were discarded by filters but were present in reference file." << endl;
			cout << "Your sample contains " << inref << " / " << ref.size() << " (" << ((double) inref / (double) ref.size()) * 100.0 << "%) junctions from the reference." << endl << endl;
		}
	}
	printFilteringResults(originalJuncs.getJunctions(),
						  filteredJuncs.getJunctions(),
						  discardedJuncs.getJunctions(),
						  string("Overall results"));
	cout << endl << "Saving junctions passing filter to disk:" << endl;
	filteredJuncs.saveAll(outputDir.string() + "/" + outputPrefix + ".pass", source + "_pass", true, this->outputExonGFF, this->outputIntronGFF);
	if (saveBad) {
		cout << "Saving junctions failing filter to disk:" << endl;
		discardedJuncs.saveAll(outputDir.string() + "/" + outputPrefix + ".fail", source + "_fail", true, this->outputExonGFF, this->outputIntronGFF);
		if (!referenceFile.empty()) {
			cout << "Saving junctions failing filters but present in reference:" << endl;
			refKeptJuncs.saveAll(outputDir.string() + "/" + outputPrefix + ".ref", source + "_ref", true, this->outputExonGFF, this->outputIntronGFF);
		}
	}
}

void portcullis::JunctionFilter::undersample(JunctionList& jl, size_t size) {
	std::mt19937 rng(12345);
	while (jl.size() > size) {
		std::uniform_int_distribution<int> gen(0, jl.size()); // uniform, unbiased
		int i = gen(rng);
		jl.erase(jl.begin() + i);
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
		for (auto & j : pass) {
			if (!j->isGenuine()) tn++;
			else fn++;
		}
		for (auto & j : fail) {
			if (j->isGenuine()) tp++;
			else fp++;
		}
	}
	else {
		for (auto & j : pass) {
			if (j->isGenuine()) tp++;
			else fp++;
		}
		for (auto & j : fail) {
			if (!j->isGenuine()) tn++;
			else fn++;
		}
	}
	return make_shared<Performance>(tp, tn, fp, fn);
}

void portcullis::JunctionFilter::forestPredict(const JunctionList& all, JunctionList& pass, JunctionList& fail, ModelFeatures& mf) {
	cout << "Creating feature vector" << endl;
	Data* testingData = mf.juncs2FeatureVectors(all);
    if (saveFeatures) {
        path feature_file = output.string() + ".features.testing";
        ofstream fout(feature_file.c_str(), std::ofstream::out);
        fout << testingData->getHeader() << endl;
        for(size_t i = 0; i < testingData->getNumRows(); i++) {
            fout << testingData->getRow(i) << endl;
        }
        fout.close();
    }

	cout << "Initialising random forest" << endl;
	shared_ptr<Forest> f = make_shared<ForestProbability>();
	vector<string> catVars;
	f->init(
		"Genuine", // Dependant variable name
		MEM_DOUBLE, // Memory mode
		testingData, // Data object
		0, // M Try (0 == use default)
		"", // Output prefix
		250, //DEFAULT_SELFTRAIN_TREES,    // Number of trees (will be overwritten when loading the model)
		1234567890, // Seed for random generator
		threads, // Number of threads
		IMP_GINI, // Importance measure
		train ? DEFAULT_MIN_NODE_SIZE_PROBABILITY : DEFAULT_MIN_NODE_SIZE_CLASSIFICATION, // Min node size
		//DEFAULT_MIN_NODE_SIZE_CLASSIFICATION,
		"", // Status var name
		true, // Prediction mode
		true, // Replace
		catVars, // Unordered categorical variable names (vector<string>)
		false, // Memory saving
		DEFAULT_SPLITRULE, // Split rule
		false, // predall
		1.0); // Sample fraction
	f->setVerboseOut(&cerr);
	// Load trees from saved model
	f->loadFromFile(modelFile.string());
	cout << "Making predictions" << endl;
	f->run(verbose);
	// Make sure score is saved back with the junction
	for (size_t i = 0; i < all.size(); i++) {
        double score = 1.0 - f->getPredictions()[i][0];
        all[i]->setScore(score);
        //cout << score << endl;
	}
	if (!genuineFile.empty() && exists(genuineFile)) {
		vector<double> thresholds;
		for (double i = 0.0; i <= 1.0; i += 0.01) {
			thresholds.push_back(i);
		}
		double max_mcc = 0.0;
		double max_f1 = 0.0;
		double best_t_mcc = 0.0;
		double best_t_f1 = 0.0;
		cout << "Threshold\t" << Performance::longHeader() << endl;
		for (auto & t : thresholds) {
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
		//threshold = best_t_mcc;
	}
	//threshold = calcGoodThreshold(f, all);
	cout << "Threshold set at " << threshold << endl;
	categorise(f, all, pass, fail, threshold);
    delete testingData;
}

void portcullis::JunctionFilter::categorise(shared_ptr<Forest> f, const JunctionList& all, JunctionList& pass, JunctionList& fail, double t) {
	for (size_t i = 0; i < all.size(); i++) {
		if ((1.0 - f->getPredictions()[i][0]) >= t) {
			pass.push_back(all[i]);
		}
		else {
			fail.push_back(all[i]);
		}
	}
}

double portcullis::JunctionFilter::calcGoodThreshold(shared_ptr<Forest> f) {
	uint32_t pos = 0;
	uint32_t neg = 0;
	for (auto & p : f->getPredictions()) {
		if ((1.0 - p[0]) >= 0.5) {
			pos++;
		}
		else {
			neg++;
		}
	}
	return (double) pos / (double) (pos + neg);
}

int portcullis::JunctionFilter::main(int argc, char *argv[]) {
	path junctionFile;
	path prepDir;
	path modelFile;
	path genuineFile;
	path filterFile;
	path referenceFile;
	path output;
	uint16_t threads;
	bool no_ml;
	bool saveBad;
    bool save_features;
	bool exongff;
	bool introngff;
	uint32_t max_length;
	string canonical;
	uint32_t mincov;
	string source;
	bool no_smote;
	bool enn;
	double threshold;
	bool verbose;
	bool help;
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	// Declare the supported options.
	po::options_description system_options("System options", w.ws_col, w.ws_col / 1.5);
	system_options.add_options()
	("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_FILTER_THREADS),
	 "The number of threads to use during testing (only applies if using forest model).")
	("verbose,v", po::bool_switch(&verbose)->default_value(false),
	 "Print extra information")
	("help", po::bool_switch(&help)->default_value(false), "Produce help message")
	;
	po::options_description output_options("Output options", w.ws_col, w.ws_col / 1.5);
	output_options.add_options()
	("output,o", po::value<path>(&output)->default_value(DEFAULT_FILTER_OUTPUT),
	 "Output prefix for files generated by this program.")
	("save_bad,b", po::bool_switch(&saveBad)->default_value(false),
	 "Saves bad junctions (i.e. junctions that fail the filter), as well as good junctions (those that pass)")
	("exon_gff", po::bool_switch(&exongff)->default_value(false),
	 "Output exon-based junctions in GFF format.")
	("intron_gff", po::bool_switch(&introngff)->default_value(false),
	 "Output intron-based junctions in GFF format.")
	("source", po::value<string>(&source)->default_value(DEFAULT_FILTER_SOURCE),
	 "The value to enter into the \"source\" field in GFF files.")
	;
	po::options_description filtering_options("Filtering options", w.ws_col, w.ws_col / 1.5);
	filtering_options.add_options()
	("filter_file,f", po::value<path>(&filterFile),
	 "If you wish to custom rule-based filter the junctions file, use this option to provide a list of the rules you wish to use.  By default we don't filter using a rule-based method, we instead filter via a self-trained random forest model.  See manual for more details.")
	("reference,r", po::value<path>(&referenceFile),
	 "Reference annotation of junctions in BED format.  Any junctions found by the junction analysis tool will be preserved if found in this reference file regardless of any other filtering criteria.  If you need to convert a reference annotation from GTF or GFF to BED format portcullis contains scripts for this.")
	("no_ml,n", po::bool_switch(&no_ml)->default_value(false),
	 "Disables machine learning filtering")
	("max_length", po::value<uint32_t>(&max_length)->default_value(0),
	 "Filter junctions longer than this value.  Default (0) is to not filter based on length.")
	("canonical", po::value<string>(&canonical)->default_value("OFF"),
	 "Keep junctions based on their splice site status.  Valid options: OFF,C,S,N. Where C = Canonical junctions (GT-AG), S = Semi-canonical junctions (AT-AC, or GC-AG), N = Non-canonical.  OFF means, keep all junctions (i.e. don't filter by canonical status).  User can separate options by a comma to keep two categories.")
	("min_cov", po::value<uint32_t>(&mincov)->default_value(1),
	 "Only keep junctions with a number of split reads greater than or equal to this number")
	("threshold", po::value<double>(&threshold)->default_value(DEFAULT_FILTER_THRESHOLD),
	 "The threshold score at which we determine a junction to be genuine or not.  Increase value towards 1.0 to increase precision, decrease towards 0.0 to increase sensitivity.  We generally find that increasing sensitivity helps when using high coverage data, or when the aligner has already performed some form of junction filtering.")
	;
	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	po::options_description hidden_options("Hidden options");
	hidden_options.add_options()
	("prep_data_dir", po::value<path>(&prepDir), "Path to directory containing prepared data.")
	("junction_file", po::value<path>(&junctionFile), "Path to the junction tab file to process.")
	("no_smote", po::bool_switch(&no_smote)->default_value(false),
	 "Use this flag to disable synthetic oversampling")
	("enn", po::bool_switch(&enn)->default_value(false),
	 "Use this flag to enable Edited Nearest Neighbour to clean decision region")
	("genuine,g", po::value<path>(&genuineFile),
	 "If you have a list of line separated boolean values in a file, indicating whether each junction in your input is genuine or not, then we can use that information here to gauge the accuracy of the predictions. This option is only useful if you have access to simulated data.")
	("model_file,m", po::value<path>(&modelFile),
	 "If you wish to use a custom random forest model to filter the junctions file, rather than self-training on the input dataset use this option to. See manual for more details.")
    ("save_features", po::bool_switch(&save_features)->default_value(false),
     "Use this flag to save features (both for training set and again for all junctions) to disk.")
	;
	// Positional option for the input bam file
	po::positional_options_description p;
	p.add("prep_data_dir", 1);
	p.add("junction_file", 1);
	// Options to display to the user
	po::options_description display_options;
	display_options.add(system_options).add(output_options).add(filtering_options);
	// Combine non-positional options
	po::options_description cmdline_options;
	cmdline_options.add(display_options).add(hidden_options);
	// Parse command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);
	// Output help information the exit if requested
	if (help || argc <= 1) {
		cout << title() << endl << endl
			 << description() << endl << endl
			 << "Usage: " << usage() << endl
			 << display_options << endl;
		return 1;
	}
	auto_cpu_timer timer(1, "\nPortcullis junction filter completed.\nTotal runtime: %ws\n\n");
	cout << "Running portcullis in junction filter mode" << endl
		 << "------------------------------------------" << endl << endl;
	// Create the prepare class
	JunctionFilter filter(prepDir, junctionFile, output);
	filter.setSaveBad(saveBad);
	filter.setSource(source);
	filter.setVerbose(verbose);
	filter.setThreads(threads);
	filter.setMaxLength(max_length);
	filter.setCanonical(canonical);
	filter.setMinCov(mincov);
	filter.setOutputExonGFF(exongff);
	filter.setOutputIntronGFF(introngff);
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
    filter.setSaveFeatures(save_features);
	filter.setReferenceFile(referenceFile);
	filter.setThreshold(threshold);
	filter.setSmote(!no_smote);
	filter.setENN(enn);
	filter.filter();
	return 0;
}

path portcullis::JunctionFilter::dataDir = ".";
