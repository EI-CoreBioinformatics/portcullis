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

#include <portcullis/bam/genome_mapper.hpp>
#include <portcullis/ml/markov_model.hpp>
#include <portcullis/junction.hpp>
using portcullis::bam::GenomeMapper;
using portcullis::ml::MarkovModel;
using portcullis::Junction;
using portcullis::JunctionPtr;
using portcullis::JunctionList;
using portcullis::SplicingScores;

namespace portcullis {
namespace ml {

typedef shared_ptr<Forest> ForestPtr;

// List of variable names
const vector<string> VAR_NAMES = {
	"Genuine",
	"rna_usrs",
	"rna_dist",
	"rna_rel",
	"rna_entropy",
	"rna_rel2raw",
	"rna_maxminanc",
	"rna_maxmmes",
	"rna_missmatch",
	"rna_intron",
	"dna_minhamm",
	"dna_coding",
	"dna_pws",
	"dna_ss"
};

struct Feature {
	string name;
	bool active;
};

class ModelFeatures {
private:
	size_t fi;
protected:
	void setRow(Data* d, size_t row, JunctionPtr j, bool labelled);

public:
	uint32_t L95;
	KmerMarkovModel exonModel;
	KmerMarkovModel intronModel;
	KmerMarkovModel donorTModel;
	KmerMarkovModel donorFModel;
	KmerMarkovModel acceptorTModel;
	KmerMarkovModel acceptorFModel;
	PosMarkovModel donorPWModel;
	PosMarkovModel acceptorPWModel;
	GenomeMapper gmap;
	vector<Feature> features;

	ModelFeatures();

	bool isCodingPotentialModelEmpty() {
		return exonModel.size() == 0 || intronModel.size() == 0;
	}

	bool isPWModelEmpty() {
		return donorPWModel.size() == 0 || acceptorPWModel.size() == 0;
	}

	void initGenomeMapper(const path& genomeFile);

	uint32_t calcIntronThreshold(const JunctionList& juncs);

	void trainCodingPotentialModel(const JunctionList& in);

	void trainSplicingModels(const JunctionList& pass, const JunctionList& fail);

	Data* juncs2FeatureVectors(const JunctionList& x);
	Data* juncs2FeatureVectors(const JunctionList& xl, const JunctionList& xu);


	ForestPtr trainInstance(const JunctionList& pos, const JunctionList& neg, string outputPrefix,
							uint16_t trees, uint16_t threads, bool probabilityMode, bool verbose, bool smote, bool enn);

	void resetActiveFeatureIndex() {
		fi = 0;
	}
	int16_t getNextActiveFeatureIndex() {
		for (int16_t i = fi + 1; i < features.size(); i++) {
			if (features[i].active) {
				fi = i;
				return i;
			}
		}
		return -1;
	}

};
}
}