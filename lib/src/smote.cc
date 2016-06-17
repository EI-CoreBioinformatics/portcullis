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

#include <portcullis/ml/smote.hpp>

portcullis::ml::Smote::Smote(uint16_t defaultK, uint16_t _smoteness, uint16_t _threads, double* _data, size_t _rows, size_t _cols) {
    
    data = _data;
    rows = _rows;
    cols = _cols;
    if (_rows < defaultK && _rows < 100)
        k = _rows;
    else
        k = defaultK;
    smoteness = _smoteness < 1 ? 1 : _smoteness;
    threads = _threads;
    verbose = false;
    
    s_rows = smoteness * rows;
    synthetic = new double[s_rows * cols];
}

void portcullis::ml::Smote::execute() {

    auto_cpu_timer timer(1, "SMOTE Time taken: %ws\n\n");

    if (verbose) {
        cout << "Starting Synthetic Minority Oversampling Technique (SMOTE)" << endl;
    }
    
    uint32_t new_index = 0;
    
    KNN knn(k, threads, data, rows, cols);
    knn.setVerbose(verbose);
    knn.execute();
    
    std::mt19937 rng(12345);
    std::uniform_int_distribution<uint16_t> igen(0, k);
    std::uniform_real_distribution<double> dgen(0, 1);
    
    
    for(size_t i = 0; i < rows; i++) {
        uint16_t N = smoteness;
        while(N > 0) {
            const vector<uint32_t> nns = knn.getNNs(i);
            uint32_t nn = nns[igen(rng)];    // Nearest neighbour row index

            for(size_t j = 0; j < cols; j++) {
                double dif = data[(nn * cols) + j] - data[(i * cols) + j];
                double gap = dgen(rng);
                synthetic[(new_index * cols) + j] = data[(i * cols) + j] + gap * dif;
            }

            new_index++;
            N--;
        }
    }    
}

void portcullis::ml::Smote::print(ostream& out) const {
    
    out << "Input:" << endl;
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            out << data[(i * cols) + j] << " ";
        }
        out << endl;
    }
    
    out << endl << "Synthetic:" << endl;
    for(size_t i = 0; i < s_rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            out << synthetic[(i * cols) + j] << " ";
        }
        out << endl;
    }
    out << endl;
}