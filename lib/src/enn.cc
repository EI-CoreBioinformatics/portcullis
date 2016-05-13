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

void portcullis::ml::ENN::doSlice(uint16_t slice, vector<bool>& result) const {

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

    // Mark for keeping or removal
    const uint16_t threshold = k / 2;
    for (uint32_t testidx = start; testidx <= end; testidx++) {

        bool pos = labels[testidx];
        uint16_t pos_count = 0;
        uint16_t neg_count = 0;
        for (size_t i = 0; i < k; i++) {
            uint32_t index = rbuf[((testidx - start) * k) + i];
            if (labels[index]) {
                pos_count++;
            } else {
                neg_count++;
            }
        }

        result[testidx] = ((pos && pos_count > threshold) || (!pos && neg_count > threshold));
    }
}

uint32_t portcullis::ml::ENN::execute(vector<bool>& results) const {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    if (verbose) {
        cout << "Performing Wilson's Edited Nearest Neighbour (ENN) ...";
        cout.flush();
    }

    results.clear();
    results.resize(rows);
    std::fill(results.begin(), results.end(), false);

    vector<thread> t(threads);
    for (uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&ENN::doSlice, this, i, std::ref(results));
    }
    for (uint16_t i = 0; i < threads; i++) {
        t[i].join();
    }

    uint32_t discard = 0;
    for (auto r : results) {
        if (!r) discard++;
    }
    uint32_t keep = results.size() - discard;

    if (verbose) {
        cout << "done." << endl << "Marked " << keep << " to be kept and " << discard << " to be discarded." << endl;
    }

    return discard;
}
