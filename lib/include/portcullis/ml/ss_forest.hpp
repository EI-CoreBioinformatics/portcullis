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

#include <memory>
using std::shared_ptr;

#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

#include <ranger/Data.h>
#include <ranger/Forest.h>

#include <portcullis/ml/model_features.hpp>
#include <portcullis/ml/markov_model.hpp>
using portcullis::ml::ModelFeatures;
using portcullis::ml::MarkovModel;
using portcullis::ml::ForestPtr;

#include <portcullis/junction.hpp>
using portcullis::bam::GenomeMapper;

using portcullis::Junction;
using portcullis::JunctionPtr;
using portcullis::JunctionList;
using portcullis::SplicingScores;

namespace portcullis {
namespace ml {
    
typedef boost::error_info<struct SSRFError,string> SSRFErrorInfo;
struct SSRFException: virtual boost::exception, virtual std::exception {};

const uint16_t REPEAT_LIMIT = 3;
    
class SemiSupervisedForest {
private:
    uint16_t trees;
    Data* labelled;
    Data* unlabelled;
    Data* all;
    ModelFeatures mf;
    ForestPtr forest;
    bool verbose;
    uint16_t threads;
    string outputPrefix;
    
public:
    SemiSupervisedForest(ModelFeatures& _mf, const JunctionList& labelled, const JunctionList& unlabelled,
            string outputPrefix, uint16_t trees, uint16_t threads, bool verbose);
    
    virtual ~SemiSupervisedForest();
    
    ForestPtr train();
    
    ForestPtr getForest() { return forest; };
    
    static Data* juncs2FeatureVectors(const JunctionList& x);
    
};
}
}