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
#include <unordered_map>
#include <vector>
using std::string;
using std::unordered_map;
using std::vector;


namespace portcullis {
namespace ml {

typedef boost::error_info<struct MMError,string> MMErrorInfo;
struct MMException: virtual boost::exception, virtual std::exception {};

/**
 * This datastructure consists of a map of kmers to a map of nucleotides with assoicated
 * probabilities.
 */
typedef unordered_map<string, unordered_map<string, double>> KMMU;
typedef unordered_map<uint32_t, unordered_map<string, double>> PMMU;

/**
 * Simple Markov chain implementation derived originally from Truesight.
 * Constructor trains the model on a set of sequences.  "getScore" returns the score
 * for a given sequence based on the pre-trained model.  
 * Automatically converts sequences to uppercase
 */
class MarkovModel {

protected:
    uint16_t order;
       
public:
    
    MarkovModel() : MarkovModel(1) {}
    
    MarkovModel(const uint32_t _order) {
        order = _order;
    }
    
    MarkovModel(const vector<string>& input, const uint32_t _order) {
        order = _order;        
        train(input);
    }
    
    void train(const vector<string>& input) {
        train(input, order);
    }
    
    virtual void train(const vector<string>& input, const uint32_t order) = 0;
    
    uint16_t getOrder() const {
        return order;
    }
    
    virtual double getScore(const string& seq) = 0;
};

class KmerMarkovModel : public portcullis::MarkovModel {
private:
    portcullis::KMMU model;
public:
    KmerMarkovModel() : MarkovModel(1) {}
    KmerMarkovModel(const uint32_t _order) : MarkovModel(_order) {};    
    KmerMarkovModel(const vector<string>& input, const uint32_t _order) : MarkovModel(input, _order) {}
    
    void train(const vector<string>& input, const uint32_t order);
    double getScore(const string& seq);
    size_t size() const {
        return model.size();
    }
};

class PosMarkovModel : public portcullis::MarkovModel {
private:
    portcullis::PMMU model;
public:
    PosMarkovModel() : MarkovModel(1) {}
    PosMarkovModel(const uint32_t _order) : MarkovModel(_order) {};    
    PosMarkovModel(const vector<string>& input, const uint32_t _order) : MarkovModel(input, _order) {}
    
    void train(const vector<string>& input, const uint32_t order);
    double getScore(const string& seq);
    size_t size() const {
        return model.size();
    }    
};

}
}

