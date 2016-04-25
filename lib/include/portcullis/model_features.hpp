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

#include <portcullis/bam/genome_mapper.hpp>
#include <portcullis/markov_model.hpp>
using portcullis::bam::GenomeMapper;
using portcullis::MarkovModel;

namespace portcullis {
    
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
};
}