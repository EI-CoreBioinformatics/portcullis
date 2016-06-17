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

#include <memory>
#include <iostream>
using std::ostream;
using std::cout;
using std::endl;
using std::make_shared;

#include <portcullis/ml/knn.hpp>

portcullis::ml::KNN::KNN(uint16_t defaultK, uint16_t _threads, double* _data, size_t _rows, size_t _cols) {
    
    data = _data;
    rows = _rows;
    cols = _cols;
    if (_rows < defaultK && _rows < 100)
        k = _rows;
    else
        k = defaultK;
    threads = _threads;
    verbose = false;
    results.resize(_rows);
}

void portcullis::ml::KNN::doSlice(uint16_t slice) {

    // Get coordinates of entries to search through in this slice
    uint32_t slice_size = rows / threads;
    uint32_t start = slice_size * slice;
    uint32_t end = (slice_size * slice) + slice_size - 1;
    if (slice == threads - 1 && end <= rows) {
        end = rows - 1; // Make sure we get the last few entries
    }
    uint32_t testcount = end - start + 1;

    vector<double> dbuf(testcount*k, std::numeric_limits<double>::max());
    vector<uint32_t> rbuf(testcount*k, 0);

    // Get feature space values for KNN for each sample in slice
    for (uint32_t baseidx = 0; baseidx < rows; baseidx++) {
        for (uint32_t testidx = start; testidx <= end; testidx++) {

            const uint32_t ri = (testidx - start) * k;

            // Get sum of squared differences (no need to do the sqrt to get 
            // the Euclidean distance... this saves about 20% runtime)
            double s = 0.0;
            for (uint16_t i = 0; i < cols; i++) {
                s += std::pow(data[baseidx * cols + i] - data[testidx * cols + i], 2.0);
            }

            // Find position to add entry...
            int i = k;
            for (; i > 0; i--)
                if (s >= dbuf[ri + i - 1])
                    break;

            // ... or skip if not if this pair is not in the current set of KNN 
            if (i >= k)
                continue;

            // Shift back entries after i
            for (int j = k - 2; j >= i; j--) {
                dbuf[ri + j + 1] = dbuf[ri + j];
                rbuf[ri + j + 1] = rbuf[ri + j];
            }

            // Record KNN and index
            dbuf[ri + i] = s;
            rbuf[ri + i] = baseidx;
        }
    }

    // Store final results
    for (uint32_t testidx = start; testidx <= end; testidx++) {
        shared_ptr<vector<uint32_t>> knn = make_shared<vector<uint32_t>>();
        
        for(size_t i = 0; i < k; i++) {
            uint32_t index = rbuf[((testidx - start) * k) + i];
            knn->push_back(index);
        }
        
        results[testidx] = knn;
    }
}

void portcullis::ml::KNN::execute() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n");

    if (verbose) {
        cout << "Performing K Nearest Neighbour (KNN) ...";
        cout.flush();
    }

    vector<thread> t(threads);
    for (uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&KNN::doSlice, this, i);
    }
    for (uint16_t i = 0; i < threads; i++) {
        t[i].join();
    }    
}

void portcullis::ml::KNN::print(ostream& out) {
    for(auto& i : results) {
        for(auto j : *i) {
            out << j << " ";
        }
        out << endl;
    }
}