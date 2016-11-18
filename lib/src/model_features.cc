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

#include <portcullis/ml/enn.hpp>
#include <portcullis/ml/smote.hpp>
using portcullis::ml::ENN;
using portcullis::ml::Smote;

#include <portcullis/junction.hpp>
using portcullis::Junction;

#include <portcullis/ml/model_features.hpp>

#include "portcullis/junction_system.hpp"



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

void portcullis::ml::ModelFeatures::setRow(Data* d, size_t row, JunctionPtr j, bool labelled) {
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
            d->set(i++, row, j->calcJunctionOverhangLogDeviation(joi), error);
        }
    }
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
        setRow(d, row, j, true);        
        row++;
    }    
    
    return d;
}

Data* portcullis::ml::ModelFeatures::juncs2FeatureVectors(const JunctionList& xl, const JunctionList& xu) {
        
    vector<string> headers;
    for(auto& f : features) {
        if (f.active) {
            headers.push_back(f.name);
        }
    }
    
    // Convert junction list info to double*
    Data* d = new DataDouble(headers, xl.size() + xu.size(), headers.size());
    
    uint32_t row = 0;
    for (const auto& j : xl) {        
        setRow(d, row, j, true);        
        row++;
    }
    
    for (const auto& j : xu) {        
        setRow(d, row, j, false);        
        row++;
    }    
    
    return d;
}



portcullis::ml::ForestPtr portcullis::ml::ModelFeatures::trainInstance(const JunctionList& pos, const JunctionList& neg,
        string outputPrefix, uint16_t trees, uint16_t threads, bool probabilityMode, bool verbose, bool smote, bool enn) {
    
    // Work out number of times to duplicate negative set
    const int N = (pos.size() / neg.size()) - 1;
    
    // Duplicate pointers to negative set
    JunctionList neg2;
    neg2.reserve(neg.size());
    neg2.insert(neg2.end(), neg.begin(), neg.end());
    
    uint32_t smote_rows = 0;
    double* smote_data = 0;
    if (N > 0 && smote) {
    
        cout << "Oversampling negative set to balance with positive set using SMOTE" << endl;
        Data* negData = juncs2FeatureVectors(neg);
        const int SC = negData->getNumCols() - 1;
    
        size_t nelements = negData->getNumRows() * SC ;
        double* nm = new double[nelements];
        for( uint32_t baseidx = 0; baseidx < negData->getNumRows(); baseidx++ ) {        
            double* r = &nm[baseidx * SC];
            for( size_t c = 1; c < negData->getNumCols(); c++) {
                r[c - 1] = negData->get(baseidx, c);
                //cout << r[c-1];
            }
            //cout << endl;
        }

        Smote smote(5, N, threads, nm, negData->getNumRows(), negData->getNumCols()-1);
        smote.execute();
        smote_rows = smote.getNbSynthRows();
        smote_data = new double[smote_rows * SC];
        double* sd = smote.getSynthetic();
        for(size_t i = 0; i < smote_rows * SC; i++) {
            smote_data[i] = sd[i];
        }
        cout << "Number of synthesized entries: " << smote.getNbSynthRows() << endl;
    }
    else if (N <= 0 && smote) {
        
        cout << "Undersampling negative set to balance with positive set" << endl;
        std::mt19937 rng(12345);
        while(neg2.size() > pos.size()) {        
            std::uniform_int_distribution<int> gen(0, neg2.size()); // uniform, unbiased
            int i = gen(rng);
            neg2.erase(neg2.begin()+i);
        }
    }
    
    
    if (verbose) cout << endl << "Combining positive, negative " << (N > 0 ? "and synthetic negative " : "") << "datasets." << endl;
    
    JunctionList training;
    training.reserve(pos.size() + neg2.size());
    training.insert(training.end(), pos.begin(), pos.end());
    training.insert(training.end(), neg2.begin(), neg2.end());

    JunctionSystem trainingSystem(training);
    trainingSystem.sort();
    JunctionList x = trainingSystem.getJunctions();
       
    
    Data* otd = juncs2FeatureVectors(x);
    
    // Create data to correct size
    Data* trainingData = N > 0 && smote ? 
        new DataDouble(
            otd->getVariableNames(), 
            x.size() + smote_rows, 
            otd->getNumCols())
        : otd;
    
    if (N > 0 && smote) {
        const int SC = trainingData->getNumCols() - 1;
        bool error = false;
        for(size_t i = 0; i < otd->getNumRows(); i++) {
            for(size_t j = 0; j < trainingData->getNumCols(); j++) {            
                trainingData->set(j, i, otd->get(i, j), error);
            }            
        }
        
        size_t k = x.size();    
        for(size_t i = 0; i < smote_rows; i++) {
            trainingData->set(0, k, 0.0, error); // Set genuine (i.e. not genuine) flag
            for(size_t j = 1; j < trainingData->getNumCols(); j++) {            
                trainingData->set(j, k, smote_data[(i * SC) + j-1], error);                
            }            
            k++;
        }
        delete[] smote_data;
    }

    vector<bool> results;
        
    if (enn) {
        if (verbose) cout << endl << "Converting training data for ENN" << endl;
        size_t elements = trainingData->getNumRows()*(trainingData->getNumCols()-1);
        double* m = new double[elements];
        for( uint32_t baseidx = 0; baseidx < trainingData->getNumRows(); baseidx++ ) {        
            double* r = &m[baseidx * (trainingData->getNumCols()-1)];
            for( size_t c = 1; c < trainingData->getNumCols(); c++) {
                r[c - 1] = trainingData->get(baseidx, c);
            }
        }    

        if (verbose) cout << "Extracting labels for ENN" << endl;
        vector<bool> labels;
        uint32_t p = 0, n = 0, o = 0;
        for(size_t i = 0; i < trainingData->getNumRows(); i++) {
            labels.push_back(trainingData->get(i, 0) == 1.0);
            if (trainingData->get(i, 0) == 1.0) {
                p++;
            }
            else if (trainingData->get(i, 0) == 0.0) {
                n++;
            }
            else {
                o++;
            }
        }
        cout << "P: " << p << "; N: " << n << "; O: " << o << endl;

        cout << endl << "Starting Wilson's Edited Nearest Neighbour (ENN) to clean decision region" << endl;
        ENN enn(3, threads, m, trainingData->getNumRows(), trainingData->getNumCols()-1, labels);
        enn.setThreshold(3);
        enn.setVerbose(true);
        uint32_t count = enn.execute(results);

        delete[] m;

        uint32_t pcount = 0, ncount = 0;
        JunctionList x2;
        for(size_t i = 0; i < trainingData->getNumRows(); i++) {
            if (trainingData->get(i, 0) == 1.0 && !results[i]) {
                pcount++;
            }
            else if (trainingData->get(i, 0) == 0.0 && !results[i]) {
                ncount++;
            }
        }

        cout << "Should discard " << pcount << " + entries and " << ncount << " - entries (Total=" << count << ")" << endl << endl;
    }
    
    uint32_t pcount = 0, ncount = 0;
    
    Data* trainingData2 = enn ? 
        new DataDouble(
            trainingData->getVariableNames(), 
            trainingData->getNumRows() - pcount,
            trainingData->getNumCols()) :
        trainingData;
    
    if (enn) {
        size_t new_index = 0;
        for(size_t i = 0; i < trainingData->getNumRows(); i++) {
            if (results[i]) { // || trainingData->get(i, 0) == 0.0) {
                bool error = false;
                for(size_t j = 0; j < trainingData->getNumCols(); j++) {            
                    trainingData2->set(j, new_index, trainingData->get(i, j), error);
                }
                new_index++;
                if (trainingData->get(i, 0) == 1) {
                    pcount++;
                }
                else {
                    ncount++;
                }
            }
        }

        delete trainingData;

        cout << "Final training set contains " << pcount << " positive entries and " << ncount << " negative entries" << endl;
    }
    
    /*path feature_file = outputPrefix + ".features";
    if (verbose) cout << "Saving feature vector to disk: " << feature_file << endl;
    
    ofstream fout(feature_file.c_str(), std::ofstream::out);    
    fout << Intron::locationOutputHeader() << "\t" << trainingData2->getHeader() << endl;
    for(size_t i = 0; i < x2.size(); i++) {        
        fout << *(x2[i]->getIntron()) << "\t" << trainingData2->getRow(i) << endl;
    }    
    fout.close();*/
         
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
        trainingData2,               // Data object
        0,                          // M Try (0 == use default)
        outputPrefix,               // Output prefix 
        trees,                      // Number of trees
        1236456789,                // Use fixed seed to avoid non-deterministic behaviour as much as possible
        threads,                    // Number of threads
        IMP_GINI,                   // Importance measure 
        probabilityMode ? DEFAULT_MIN_NODE_SIZE_PROBABILITY : DEFAULT_MIN_NODE_SIZE_CLASSIFICATION,  // Min node size
        "",                         // Status var name 
        false,                      // Prediction mode
        false,                       // Replace 
        catVars,                    // Unordered categorical variable names (vector<string>)
        false,                      // Memory saving
        AUC, //DEFAULT_SPLITRULE,          // Split rule
        false,                      // predall
        1.0);                       // Sample fraction
            
    if (verbose) cout << "Training" << endl;
    f->setVerboseOut(&cerr);
    f->run(verbose);
    cout << "OOBE: " << f->getOverallPredictionError() << endl;
    
    delete trainingData2;
    
    return f;
}