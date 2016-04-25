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

#include <iostream>
using std::cout;
using std::endl;

#include <ranger/DataDouble.h>

#include <portcullis/junction.hpp>
using portcullis::Junction;

#include <portcullis/model_features.hpp>



void portcullis::ModelFeatures::initGenomeMapper(const path& genomeFile) {

    // Initialise
    gmap.setGenomeFile(genomeFile);

    // Load the fasta index
    gmap.loadFastaIndex();
}

uint32_t portcullis::ModelFeatures::calcIntronThreshold(const JunctionList& juncs) {

    vector<uint32_t> intron_sizes;
    for(auto& j : juncs) {
        intron_sizes.push_back(j->getIntronSize());
    }

    std::sort(intron_sizes.begin(), intron_sizes.end());

    L95 = intron_sizes[intron_sizes.size() * 0.95];

    return L95;
}

void portcullis::ModelFeatures::trainCodingPotentialModel(const JunctionList& in) {

    vector<string> exons;
    vector<string> introns;
    for(auto& j : in) {

        int len = 0;

        string left_exon = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start - 80, j->getIntron()->start, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            left_exon = SeqUtils::reverseComplement(left_exon);
        }        
        exons.push_back(left_exon);

        string left_intron = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->start, j->getIntron()->start+81, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            left_intron = SeqUtils::reverseComplement(left_intron);
        }        
        introns.push_back(left_intron);

        string right_intron = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->end-80, j->getIntron()->end+1, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            right_intron = SeqUtils::reverseComplement(right_intron);
        }        
        introns.push_back(right_intron);


        string right_exon = gmap.fetchBases(j->getIntron()->ref.name.c_str(), j->getIntron()->end + 1, j->getIntron()->end + 80, &len);
        if (j->getConsensusStrand() == Strand::NEGATIVE) {
            right_exon = SeqUtils::reverseComplement(right_exon);
        }        
        exons.push_back(right_exon);
    }

    exonModel.train(exons, 5);
    intronModel.train(introns, 5);
}

Data* portcullis::ModelFeatures::juncs2FeatureVectors(const JunctionList& x) {
    
    vector<string> headers;
    headers.reserve( VAR_NAMES.size() + JO_NAMES.size() );
    headers.insert( headers.end(), VAR_NAMES.begin(), VAR_NAMES.end() );
    headers.insert( headers.end(), JO_NAMES.begin()+14, JO_NAMES.end() );
    
    // Convert junction list info to double*
    double* d = new double[x.size() * headers.size()];
    
    uint32_t row = 0;
    for (const auto& j : x) {        
        
        //d[0 * x.size() + row] = j->getNbJunctionAlignments();
        //d[0 * x.size() + row] = j->getNbDistinctAlignments();
        d[0 * x.size() + row] = j->getNbReliableAlignments();
        //d[3 * x.size() + row] = j->getMaxMinAnchor();
        //d[4 * x.size() + row] = j->getDiffAnchor();
        //d[5 * x.size() + row] = j->getNbDistinctAnchors();
        d[1 * x.size() + row] = j->getEntropy();
        d[2 * x.size() + row] = j->getMaxMMES();
        d[3 * x.size() + row] = std::min(j->getHammingDistance5p(), j->getHammingDistance3p());
        //d[3 * x.size() + row] = ;
        d[4 * x.size() + row] = j->getReliable2RawRatio();
        //d[5 * x.size() + row] = j->getMeanMismatches();
        d[5 * x.size() + row] = L95 == 0 ? 0.0 : j->calcIntronScore(L95);     // Intron score
        d[6 * x.size() + row] = isCodingPotentialModelEmpty() ? 0.0 : j->calcCodingPotential(gmap, exonModel, intronModel);     // Coding potential score
        d[7 * x.size() + row] = j->isGenuine();
        // Junction overhang values at each position are first converted into deviation from expected distributions       
        double half_read_length = (double)j->getMeanQueryLength() / 2.0;    
        for(size_t i = 14; i < JO_NAMES.size(); i++) {
            double Ni = j->getJunctionOverhangs(i);                 // Actual count at this position
            if (Ni == 0.0) Ni = 0.000000000001;                     // Ensure some value > 0 here otherwise we get -infinity later.
            double Pi = 1.0 - ((double)i / half_read_length);       // Likely scale at this position
            double Ei = (double)j->getNbJunctionAlignments() * Pi;  // Expected count at this position
            double Xi = abs(log2(Ni / Ei));                         // Log deviation
            
            d[(i-14 + 8) * x.size() + row] = Xi;
        }
        
        row++;
    }
    
    Data* data = new DataDouble(d, headers, x.size(), headers.size());
    data->setExternalData(false);      // This causes 'd' to be deleted when 'data' is deleted
    return data;
}