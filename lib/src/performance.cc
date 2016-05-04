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
using std::cout;
using std::endl;

#include <portcullis/performance.hpp>


string portcullis::ml::Performance::toShortString() const {
    vector<string> parts;
    parts.push_back(std::to_string(tp));
    parts.push_back(std::to_string(tn)); 
    parts.push_back(std::to_string(fp));
    parts.push_back(std::to_string(fn));
    parts.push_back(portcullis::Performance::to_2dp_string(getRecall()));
    parts.push_back(portcullis::Performance::to_2dp_string(getPrecision()));
    parts.push_back(portcullis::Performance::to_2dp_string(getF1Score()));
    return boost::algorithm::join(parts, "\t");
}

string portcullis::ml::Performance::toLongString() const {
    vector<string> parts;
    parts.push_back(std::to_string(tp));
    parts.push_back(std::to_string(tn)); 
    parts.push_back(std::to_string(fp));
    parts.push_back(std::to_string(fn));
    parts.push_back(portcullis::Performance::to_2dp_string(getPrevalence()));        
    parts.push_back(portcullis::Performance::to_2dp_string(getBias()));
    parts.push_back(portcullis::Performance::to_2dp_string(getSensitivity()));
    parts.push_back(portcullis::Performance::to_2dp_string(getSpecificity()));
    parts.push_back(portcullis::Performance::to_2dp_string(getPrecision()));
    parts.push_back(portcullis::Performance::to_2dp_string(getNPV()));
    parts.push_back(portcullis::Performance::to_2dp_string(getF1Score()));
    parts.push_back(portcullis::Performance::to_2dp_string(getAccuracy()));
    parts.push_back(portcullis::Performance::to_2dp_string(getInformedness()));
    parts.push_back(portcullis::Performance::to_2dp_string(getMarkedness()));
    parts.push_back(portcullis::Performance::to_2dp_string(getMCC()));        
    return boost::algorithm::join(parts, "\t");
}

    
void portcullis::ml::Performance::loadGenuine(path& genuineFile, vector<bool>& results) {

    // Load reference data    
    std::ifstream refs(genuineFile.string());
    string line;
    uint32_t lineNb = 0;
    while (std::getline(refs, line)) {
        std::istringstream iss(line);
        bool res;
        iss >> res;
        results.push_back(res);
    }
    refs.close();
}

void portcullis::ml::PerformanceList::outputMeanPerformance(std::ostream& resout) {
    
    vector<double> prevs;
    vector<double> biases;
    vector<double> recs;
    vector<double> prcs;
    vector<double> spcs;
    vector<double> f1s;
    vector<double> infs;
    vector<double> mrks;
    vector<double> accs;
    vector<double> mccs;
        
    for(auto& p : this->scores) {
        prevs.push_back(p->getPrevalence());
        biases.push_back(p->getBias());
        recs.push_back(p->getRecall());
        prcs.push_back(p->getPrecision());
        spcs.push_back(p->getSpecificity());
        f1s.push_back(p->getF1Score());
        infs.push_back(p->getInformedness());
        mrks.push_back(p->getMarkedness());
        accs.push_back(p->getAccuracy());
        mccs.push_back(p->getMCC());
    }
    
    outputMeanScore(prevs, "prevalence", resout);
    outputMeanScore(biases, "bias", resout);
    outputMeanScore(recs, "recall", resout);
    outputMeanScore(prcs, "precision", resout);
    outputMeanScore(f1s, "F1", resout);
    outputMeanScore(spcs, "specificity", resout);
    outputMeanScore(accs, "accuracy", resout);
    outputMeanScore(infs, "informedness", resout);
    outputMeanScore(mrks, "markededness", resout);
    outputMeanScore(mccs, "MCC", resout);
}

void portcullis::ml::PerformanceList::outputMeanScore(const vector<double>& scores, const string& score_type, std::ostream& resout) {
    double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
    double mean = sum / scores.size();
    double sq_sum = std::inner_product(scores.begin(), scores.end(), scores.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / scores.size() - mean * mean);
    
    stringstream msg;
    msg << "Mean " << std::left << std::setw(13) << score_type << ": " << std::fixed << std::setprecision(2) << mean << "% (+/- " << stdev << "%)" << endl;

    cout << msg.str();
    resout << msg.str();
}

