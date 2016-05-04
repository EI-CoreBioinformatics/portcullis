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

#include <string>
#include <algorithm>
#include <math.h>
#include <iterator>
#include <memory>
#include <vector>
using std::string;
using std::stringstream;
using std::vector;
using std::random_shuffle;
using std::shared_ptr;
using std::unique_ptr;

#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

#include <ranger/Forest.h>

#include <portcullis/ml/model_features.hpp>
using portcullis::ml::ModelFeatures;
using portcullis::ml::ForestPtr;

namespace portcullis {
    
typedef boost::error_info<struct TrainError,string> TrainErrorInfo;
struct TrainException: virtual boost::exception, virtual std::exception { };

const string DEFAULT_TRAIN_OUTPUT = "portcullis_train/portcullis";
const uint16_t DEFAULT_TRAIN_FOLDS = 5;
const uint16_t DEFAULT_TRAIN_TREES = 100;
const uint16_t DEFAULT_TRAIN_THREADS = 1;
const double DEFAULT_TRAIN_FRACTION = 1.0;
const int DEFAULT_SEED = 1234567;       // To avoid non-deterministic behaviour


class Train {
    
    
private:
    
    path junctionFile;
    path refFile;
    path outputPrefix;
    uint16_t folds;
    uint16_t trees;
    uint16_t threads;
    double fraction;
    bool regressionMode;
    bool verbose;    
    
public:

    
    
    Train() {}
    
    Train(const path& _junctionFile, const path& _refFile);
    
    virtual ~Train() {
    }
    
    uint16_t getFolds() const {
        return folds;
    }

    void setFolds(uint16_t folds) {
        this->folds = folds;
    }

    uint16_t getTrees() const {
        return trees;
    }

    void setTrees(uint16_t trees) {
        this->trees = trees;
    }
    
    uint16_t getThreads() const {
        return threads;
    }

    void setThreads(uint16_t threads) {
        this->threads = threads;
    }
    
    path getRefFile() const {
        return refFile;
    }

    void setRefFile(path refFile) {
        this->refFile = refFile;
    }

    path getJunctionFile() const {
        return junctionFile;
    }

    void setJunctionFile(path junctionFile) {
        this->junctionFile = junctionFile;
    }

    path getOutputPrefix() const {
        return outputPrefix;
    }

    void setOutputPrefix(path outputPrefix) {
        this->outputPrefix = outputPrefix;
    }

    double getFraction() const {
        return fraction;
    }

    void setFraction(double fraction) {
        this->fraction = fraction;
    }

    bool isRegressionMode() const {
        return regressionMode;
    }

    void setRegressionMode(bool regressionMode) {
        this->regressionMode = regressionMode;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }

    
    
    /**
     * Run a supervised training algorithm on the input data using k fold cross validation
     */
    void train();
    
    
    static string helpMessage() {
        return string("\nPortcullis Training Mode Help.\n\n") +
                      "This is mode is intended to train a random forest model for later use\n" +
                      "by the junction filtering tool.  The current implementation uses \"ranger\"\n" +
                      "for creating the decision trees and random forests.  We also provide k-fold\n" +
                      "cross validation support to help reduce overfitting.\n\n" +
                      "Usage: portcullis train [options] --reference=<labels_file> <junction_file>\n\n" +
                      "Allowed options";
    }

    static int main(int argc, char *argv[]);
    
    
protected:
    
    
    void testInstance(ForestPtr f, const JunctionList& y);
    
    void getRandomSubset(const JunctionList& in, JunctionList& out);    
};
}

