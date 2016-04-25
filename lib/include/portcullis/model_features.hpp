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

#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

#include <ranger/Data.h>

#include <portcullis/bam/genome_mapper.hpp>
#include <portcullis/markov_model.hpp>
#include <portcullis/junction.hpp>
using portcullis::bam::GenomeMapper;
using portcullis::MarkovModel;
using portcullis::JunctionList;

namespace portcullis {

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
            "CodingPotential",
            "Genuine" };
    
    
class ModelFeatures {
public:
    uint32_t L95;
    MarkovModel exonModel;
    MarkovModel intronModel;
    GenomeMapper gmap;
    
    ModelFeatures() : L95(0) {
    }
    
    bool isCodingPotentialModelEmpty() {
        return exonModel.size() == 0 || intronModel.size() == 0;
    }
    
    void initGenomeMapper(const path& genomeFile);
    
    uint32_t calcIntronThreshold(const JunctionList& juncs);
    
    void trainCodingPotentialModel(const JunctionList& in);
    
    Data* juncs2FeatureVectors(const JunctionList& x);
};
}