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

#include <portcullis/ml/icote.hpp>

portcullis::ml::Icote::Icote(uint16_t defaultK, uint16_t _duplications, uint16_t _threads, double* _data, size_t _rows, size_t _cols) {
    
    data = _data;
    rows = _rows;
    cols = _cols;
    if (_rows < defaultK && _rows < 100)
        k = _rows;
    else
        k = defaultK;
    duplications = _duplications < 1 ? 1 : _duplications;
    threads = _threads;
    verbose = false;
    
    s_rows = duplications * rows;
    
}

void portcullis::ml::Smote::execute() {

    auto_cpu_timer timer(1, "ICOTE Time taken: %ws\n\n");

    if (verbose) {
        cout << "Starting Immume Centroids Oversampling Technique (ICOTE)" << endl;
    }
    
    // Attribute selection
    // To reduce computational cost remove constant attributes
    if (verbose) cout << "Selecting non-constant features" << endl;
    
    vector<bool> is_const(cols);
    vector<double> const_val(cols);
    uint16_t nb_const_col = 0;
    for(size_t k = 0; k < cols; k++) {
        is_const[k] = this->isConstant(k);
        if (is_const[k]) {
            nb_const_col++;
            const_val[k] = data[k];
        }
    }
    uint16_t sel_col = cols - nb_const_col;
    
    if (verbose) cout << "Ignoring " << const_col << " columns as values are constant across all samples" << endl;
    
    double* ag = new double[rows * sel_col];
    for(size_t i = 0; i < rows; i++) {
        size_t new_col = 0;
        for(size_t j = 0; j < cols; j++) {
            if (!is_const[j]) {
                ag[(i * sel_col) + new_col++] = data[(i*cols) + j];
            }
        }
    }
    
    if (verbose) cout << "Created non-constant feature matrix" << endl;
    
    
    // Normalisation
    if (verbose) cout << "Normalising feature matrix" << endl;
    vector<double> min_col(sel_col);
    vector<double> max_col(sel_col);
    for(size_t j = 0; j < sel_col; j++) {
        double min_v = std::numeric_limits<double>::max();
        double max_v = std::numeric_limits<double>::min();
        for(size_t i = 0; i < rows; i++) {
            double x = ag[(i * sel_col) + j];
            min_v = std::min(min_v, x);
            max_v = std::max(max_v, x);
        }
        min_col[j] = min_v;
        max_col[j] = max_v;
    }
    
    for(size_t j = 0; j < sel_col; j++) {
        double min_x = min_col[j];
        double max_x = max_col[j];
        for(size_t i = 0; i < rows; i++) {
            double x = ag[(i * sel_col) + j];
            ag[(i * sel_col) + j] = (x - min_x) / (max_x - min_x);
        }
    }
    if (verbose) cout << "Normalised feature matrix" << endl;
    
    if (verbose) cout << "Generating Random Antibodies" << endl;
    vector<double> ab(k * sel_col);
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dgen(0.0, 1.0);
    
    for(size_t i = 0; i < k; i++) {
        for(size_t j = 0; j < sel_col; j++) {
            ab[(i*sel_col) + j] = dgen(rng);
        }
    }
    if (verbose) cout << "Random Antibodies generated" << endl;
    
    uint16_t N = 0;
    while (N < 10) {
        vector<double> M(duplications * rows * sel_col);
        
        for(size_t i = 0; i < rows; i++) {
            vector<double> dist(rows * duplications * sel_col);
            for(size_t j = 0; j < duplications; j++) {
                double s = 0.0;
                for (size_t c = 0; c < sel_cols; c++) {
                    s += std::pow(ag[i * sel_cols + c] - ab[i * sel_cols + c], 2.0);
                }
                dist[(i * sel_cols) + j] = std::sqrt(s);
            }
        }
    
        if (verbose) cout << "Mutating antibodies" << endl;
    
        
        double s = 0.0;
        for (size_t j = 0; j < sel_cols; j++) {
            s += std::pow(antibodies[i * sel_cols + j] - attr_sel[i * sel_cols + j], 2.0);
        }
        double ak = (1.0 / sel_cols) * std::sqrt(ak);
    }
    
    if (verbose) cout << "Antibodies mutated" << endl;
    
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