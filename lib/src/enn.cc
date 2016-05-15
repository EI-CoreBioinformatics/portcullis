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

#include <portcullis/ml/knn.hpp>
using portcullis::ml::KNN;

#include <portcullis/ml/enn.hpp>

portcullis::ml::ENN::ENN(uint16_t defaultK, uint16_t _threads, double* _data, size_t _rows, size_t _cols, vector<bool>& _labels) {
    if (_rows != _labels.size()) {
        BOOST_THROW_EXCEPTION(ENNException() << ENNErrorInfo(string(
                "The supplied number of rows does not match the number of labels")));
    }

    data = _data;
    rows = _rows;
    cols = _cols;
    labels = _labels;
    if (_rows < defaultK && _rows < 100)
        k = _rows;
    else
        k = defaultK;
    threads = _threads;
    verbose = false;
}

uint32_t portcullis::ml::ENN::execute(vector<bool>& results) const {

    auto_cpu_timer timer(1, "ENN Time taken: %ws\n\n");

    if (verbose) {
        cout << "Starting Wilson's Edited Nearest Neighbour (ENN)" << endl;
    }

    KNN knn(k, threads, data, rows, cols);
    knn.setVerbose(verbose);
    knn.execute();    
    
    if (verbose) {
        cout << "Finding outliers" << endl;
    }
    
    results.clear();
    results.resize(rows);
    uint32_t discard = 0;
    
    const uint16_t threshold = this->getThreshold();
    
    for(size_t i = 0; i < rows; i++) {
        uint16_t pos_count = 0;
        uint16_t neg_count = 0;
        bool pos = labels[i];
        const vector<uint32_t> nn = knn.getNNs(i);
        for (size_t j = 0; j < k; j++) {
            uint32_t index = nn[j];
            if (labels[index]) {
                pos_count++;
            } else {
                neg_count++;
            }
        }
        results[i] = ((pos && pos_count > threshold) || (!pos && neg_count > threshold));    
        if (!results[i]) discard++;
    }

    uint32_t keep = results.size() - discard;

    if (verbose) {
        cout << "Marked " << keep << " to be kept and " << discard << " to be discarded." << endl;
    }

    return discard;
}
