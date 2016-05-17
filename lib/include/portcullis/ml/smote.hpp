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

#include <boost/exception/all.hpp>

namespace portcullis {
namespace ml {
    
typedef boost::error_info<struct SmoteError,string> SmoteErrorInfo;
struct SmoteException: virtual boost::exception, virtual std::exception { };


/**
 * Object to perform classification on balanced ensembled selected from random sampling.
 * 
 * It is based on the idea presented in the paper "Exploratory Undersampling 
 * Class-Imbalance Learning" by Liu et al.    
 */
class Smote {
private:
    uint16_t k;
    uint16_t smoteness;
    uint16_t threads;
    bool verbose;
    
    double* data;
    size_t rows;
    size_t cols;
     
    double* synthetic;
    size_t s_rows;
        
public:

    Smote(uint16_t defaultK, uint16_t _smoteness, uint16_t _threads, double* _data, size_t _rows, size_t _cols);
    
    ~Smote() {
        delete[] synthetic;
    }
    
    uint16_t getK() const {
        return k;
    }

    void setK(uint16_t k) {
        this->k = k;
    }

    uint16_t getThreads() const {
        return threads;
    }

    void setThreads(uint16_t threads) {
        this->threads = threads;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }
    
    double* getSynthetic() {
        return synthetic;
    }
    
    size_t getNbSynthRows() const {
        return s_rows;
    }
    
    double getSynth(size_t row, size_t col) const {
        return synthetic[(row * cols) + col];
    }
    
    void execute();
    
    
};
}
}