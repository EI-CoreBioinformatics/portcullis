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

#include <fstream>
#include <iostream>
#include <memory>
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::make_shared;

#include <ranger/DataDouble.h>
#include <ranger/ForestProbability.h>
#include <ranger/ForestClassification.h>

#include <portcullis/junction.hpp>
using portcullis::Junction;

#include <portcullis/ml/model_features.hpp>



void portcullis::ml::ModelFeatures::initGenomeMapper(const path& genomeFile) {

    // Initialise
    gmap.setGenomeFile(genomeFile);

    // Load the fasta index
    gmap.loadFastaIndex();
}

uint32_t portcullis::ml::ModelFeatures::calcIntronThreshold(const JunctionList& juncs) {

    vector<uint32_t> intron_sizes;
    for(auto& j : juncs) {
        intron_sizes.push_back(j->getIntronSize());
    }

    std::sort(intron_sizes.begin(), intron_sizes.end());

    L95 = intron_sizes[intron_sizes.size() * 0.95];

    return L95;
}

void portcullis::ml::ModelFeatures::trainCodingPotentialModel(const JunctionList& in) {

    vector<string> exons;
    vector<string> introns;
    for(auto& j : in) {

        int len = 0;

        string left_exon = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start - 202, j->getIntron()->start - 2, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            left_exon = SeqUtils::reverseComplement(left_exon);
        }        
        exons.push_back(left_exon);

        /*string left_intron = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start, j->getIntron()->start+80, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            left_intron = SeqUtils::reverseComplement(left_intron);
        }        
        introns.push_back(left_intron);

        string right_intron = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->end-80, j->getIntron()->end, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            right_intron = SeqUtils::reverseComplement(right_intron);
        }        
        introns.push_back(right_intron);*/
        
        string intron = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start, j->getIntron()->end, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            intron = SeqUtils::reverseComplement(intron);
        }        
        introns.push_back(intron);


        string right_exon = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->end + 1, j->getIntron()->end + 201, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            right_exon = SeqUtils::reverseComplement(right_exon);
        }        
        exons.push_back(right_exon);
    }

    exonModel.train(exons, 5);
    intronModel.train(introns, 5);
}

void portcullis::ml::ModelFeatures::trainSplicingModels(const JunctionList& pass, const JunctionList& fail) {

    vector<string> donors;
    vector<string> acceptors;
    for(auto& j : pass) {

        int len = 0;

        string left = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start - 3, j->getIntron()->start + 20, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            left = SeqUtils::reverseComplement(left);
        }        
   
        string right = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->end - 20, j->getIntron()->end + 2, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            right = SeqUtils::reverseComplement(right);
        }        
        
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            donors.push_back(right);
            acceptors.push_back(left);
        }
        else {
            donors.push_back(left);
            acceptors.push_back(right);            
        }
        
    }

    donorPWModel.train(donors, 1);
    acceptorPWModel.train(acceptors, 1);
    donorTModel.train(donors, 5);
    acceptorTModel.train(acceptors, 5);
    
    donors.clear();
    acceptors.clear();
    
    for(auto& j : fail) {

        int len = 0;

        string left = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start - 3, j->getIntron()->start + 20, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            left = SeqUtils::reverseComplement(left);
        }        
   
        string right = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->end - 20, j->getIntron()->end + 2, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            right = SeqUtils::reverseComplement(right);
        }        
        
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            donors.push_back(right);
            acceptors.push_back(left);
        }
        else {
            donors.push_back(left);
            acceptors.push_back(right);
        }
    }

    donorFModel.train(donors, 5);
    acceptorFModel.train(acceptors, 5);
}

Data* portcullis::ml::ModelFeatures::juncs2FeatureVectors(const JunctionList& x) {
        
    vector<string> headers;
    for(auto& f : features) {
        if (f.active) {
            headers.push_back(f.name);
        }
    }
    
    // Convert junction list info to double*
    Data* d = new DataDouble(headers, x.size(), headers.size());
    
    uint32_t row = 0;
    for (const auto& j : x) {
        SplicingScores ss = j->calcSplicingScores(gmap, donorTModel, donorFModel, acceptorTModel, acceptorFModel, 
                donorPWModel, acceptorPWModel);
        
        bool error = false;
        d->set(0, row, j->isGenuine(), error);
        
        uint16_t i = 1;
        
        if (features[1].active) {
            d->set(i++, row, j->getNbUniquelySplicedReads(), error);
        }
        if (features[2].active) {
            d->set(i++, row, j->getNbDistinctAlignments(), error);
        }
        if (features[3].active) {
            d->set(i++, row, j->getNbReliableAlignments(), error);
        }
        if (features[4].active) {
            d->set(i++, row, j->getEntropy(), error);
        }
        if (features[5].active) {
            d->set(i++, row, j->getReliable2RawRatio(), error);
        }        
        if (features[6].active) {
            d->set(i++, row, j->getMaxMinAnchor(), error);
        }        
        if (features[7].active) {
            d->set(i++, row, j->getMaxMMES(), error);
        }
        if (features[8].active) {
            d->set(i++, row, j->getMeanMismatches(), error);
        }
        if (features[9].active) {
            d->set(i++, row, L95 == 0 ? 0.0 : j->calcIntronScore(L95), error);
        }
        if (features[10].active) {
            d->set(i++, row, std::min(j->getHammingDistance5p(), j->getHammingDistance3p()), error);
        }
        if (features[11].active) {
            d->set(i++, row, isCodingPotentialModelEmpty() ? 0.0 : j->calcCodingPotential(gmap, exonModel, intronModel), error);
        }
        if (features[12].active) {
            d->set(i++, row, isPWModelEmpty() ? 0.0 : ss.positionWeighting, error);
        }
        if (features[13].active) {
            d->set(i++, row, isPWModelEmpty() ? 0.0 : ss.splicingSignal, error);
        }
                
        //Junction overhang values at each position are first converted into deviation from expected distributions       
        for(size_t joi = 0; joi < JO_NAMES.size(); joi++) {
            if (features[joi + 14].active) {
                d->set(i++, row, j->getJunctionOverhangLogDeviation(joi), error);
            }
        }
        
        row++;
    }
    
    return d;
}



portcullis::ml::ForestPtr portcullis::ml::ModelFeatures::trainInstance(const JunctionList& x, 
        string outputPrefix, uint16_t trees, uint16_t threads, bool probabilityMode, bool verbose) {
    
    if (verbose) cout << "Creating feature vector" << endl;
    Data* trainingData = juncs2FeatureVectors(x);
    
    path feature_file = outputPrefix + ".features";
    if (verbose) cout << "Saving feature vector to disk: " << feature_file << endl;
    
    ofstream fout(feature_file.c_str(), std::ofstream::out);
    
    fout << Intron::locationOutputHeader() << "\t" << trainingData->getHeader() << endl;
    for(size_t i = 0; i < x.size(); i++) {        
        fout << *(x[i]->getIntron()) << "\t" << trainingData->getRow(i) << endl;
    }
    
    fout.close();
     
    if (verbose) cout << "Initialising random forest" << endl;
    ForestPtr f = nullptr;
    if (probabilityMode) {
        f = make_shared<ForestProbability>();
    }
    else {
        f = make_shared<ForestClassification>();
    }
    
    vector<string> catVars;
    
    f->init(
        "Genuine",                  // Dependant variable name
        MEM_DOUBLE,                 // Memory mode
        trainingData,               // Data object
        0,                          // M Try (0 == use default)
        outputPrefix,               // Output prefix 
        250, //trees,                      // Number of trees
        1236456789,                // Use fixed seed to avoid non-deterministic behaviour as much as possible
        threads,                    // Number of threads
        IMP_GINI,                   // Importance measure 
        probabilityMode ? DEFAULT_MIN_NODE_SIZE_PROBABILITY : DEFAULT_MIN_NODE_SIZE_CLASSIFICATION,  // Min node size
        "",                         // Status var name 
        false,                      // Prediction mode
        true,                       // Replace 
        catVars,                    // Unordered categorical variable names (vector<string>)
        false,                      // Memory saving
        DEFAULT_SPLITRULE,          // Split rule
        false,                      // predall
        1.0);                       // Sample fraction
            
    if (verbose) cout << "Training" << endl;
    f->setVerboseOut(&cerr);
    f->run(verbose);
    
    return f;
}