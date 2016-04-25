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

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>

#include <portcullis/markov_model.hpp>

void portcullis::MarkovModel::train(const vector<string>& input, const uint32_t _order) {
    order = _order;
    MMU temp;
    for(auto& seq : input) {
        string sequp = boost::to_upper_copy(seq);
        for(uint16_t i = order; i < sequp.size(); i++) {
            temp[sequp.substr(i-order, order)][sequp.substr(i, 1)]++;
        }
    }
    model.clear();
    for(auto& i : temp) {
        double sum = 0;
        for(auto& j : i.second) {
            sum += j.second;
        }
        for(auto& j : i.second) {
            //cout << i.first << " " << j.first << " " << j.second << endl;
            model[i.first][j.first] = j.second / sum;
        }
    }
}
    
    
double portcullis::MarkovModel::getScore(const string& seq) {
    string sequp = boost::to_upper_copy(seq); 
    double score = 1.0;
    for(uint16_t i = order; i < sequp.size(); i++){
        score *= model[sequp.substr(i-order, order)][sequp.substr(i, 1)];
    }
    if(score == 0.0) {
        return -100.0;
    }
    return log(score);
}
