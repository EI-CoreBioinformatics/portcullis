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
#include <portcullis/performance.hpp>


namespace portcullis {
    
typedef boost::error_info<struct TrainError,string> TrainErrorInfo;
struct TrainException: virtual boost::exception, virtual std::exception { };

const string DEFAULT_TRAIN_OUTPUT = "portcullis_train/portcullis";
const uint16_t DEFAULT_TRAIN_FOLDS = 5;
const uint16_t DEFAULT_TRAIN_TREES = 100;
const uint16_t DEFAULT_TRAIN_THREADS = 1;
const double DEFAULT_TRAIN_FRACTION = 1.0;
const int DEFAULT_SEED = 1234567;       // To avoid non-deterministic behaviour

// List of variable names
const vector<string> VAR_NAMES = { 
            //"M2-nb-reads", 
            //"M3-nb_dist_aln", 
            "nb_rel_aln", 
            //"M8-max_min_anc", 
            //"M9-dif_anc", 
            //"M10-dist_anc", 
            "entropy", 
            "maxmmes", 
            "min_hamming_score", 
            //"M14-hamming3p",
            "rel2raw_ratio",
            //"mean_mismatches",
            "IntronScore",
            "Genuine" };

// Derived from https://sureshamrita.wordpress.com/2011/08/24/c-implementation-of-k-fold-cross-validation/
template<class In>
class KFold {
public:
    KFold(int k, In _beg, In _end) :
        beg(_beg), end(_end), K(k) {
        if (K <= 0)
            BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "The supplied value of K is =") + lexical_cast<string>(K) + 
                        ". One cannot create " + lexical_cast<string>(K) + "no of folds"));

        //create the vector of integers
        int foldNo = 0;
        for (In i = beg; i != end; i++) {
            whichFoldToGo.push_back(++foldNo);
            if (foldNo == K)
                foldNo = 0;
        }
        if (!K)
            BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "With this value of k (=") + lexical_cast<string>(K) + 
                        ")Equal division of the data is not possible"));

        random_shuffle(whichFoldToGo.begin(), whichFoldToGo.end());
    }
        
    template<class Out>
    void getFold(int foldNo, Out training, Out testing) {
        int k = 0;
        In i = beg;
        while (i != end) {
            if (whichFoldToGo[k++] == foldNo) {
                *testing++ = *i++;
            } else
                *training++ = *i++;
        }
    }
    
private:
    In beg;
    In end;
    int K; //how many folds in this
    vector<int> whichFoldToGo; 
};


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
    
    static Data* juncs2FeatureVectors(const JunctionList& x) {
        return juncs2FeatureVectors(x, 0);
    }
    static Data* juncs2FeatureVectors(const JunctionList& x, const uint32_t L95);
    
    static shared_ptr<Forest> trainInstance(const JunctionList& x, string outputPrefix, 
            uint16_t trees, uint16_t threads, bool regressionMode, bool verbose) {
        return trainInstance(x, outputPrefix, trees, threads, regressionMode, verbose, 0);
    }
    static shared_ptr<Forest> trainInstance(const JunctionList& x, string outputPrefix, 
            uint16_t trees, uint16_t threads, bool regressionMode, bool verbose, const uint32_t L95);
    
protected:
    
    
    void testInstance(shared_ptr<Forest> f, const JunctionList& y) {
        return testInstance(f, y, 0);
    }
    
    void testInstance(shared_ptr<Forest> f, const JunctionList& y, const uint32_t L95);
    
    void getRandomSubset(const JunctionList& in, JunctionList& out);
    
    void outputMeanPerformance(const vector<unique_ptr<Performance>>& scores, std::ofstream& resout);
    
    void outputMeanScore(const vector<double>& scores, const string& score_type, std::ofstream& resout);
};
}

