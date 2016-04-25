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