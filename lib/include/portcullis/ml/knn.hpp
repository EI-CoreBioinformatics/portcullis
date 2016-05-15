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

#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>
using std::ostream;
using std::string;
using std::thread;
using std::vector;
using std::shared_ptr;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

#include <ranger/Data.h>

namespace portcullis {
namespace ml {

typedef boost::error_info<struct KNNError,string> KNNErrorInfo;
struct KNNException: virtual boost::exception, virtual std::exception { };
    
/**
 * An parallel implementation of K Nearest Neighbour.
 * Logic originally derived from OpenCV
 */
class KNN {
protected:
    uint16_t k;
    uint16_t threads;
    bool verbose;
    
    double* data;
    size_t rows;
    size_t cols;
    
    vector<shared_ptr<vector<uint32_t>>> results;
    
    void doSlice( uint16_t slice );
    
public:    
    
    KNN(uint16_t defaultK, uint16_t _threads, double* _data, size_t _rows, size_t _cols);
    
    int getK() const {
        return k;
    }

    void setK(int k) {
        this->k = k;
    }

    int getThreads() const {
        return threads;
    }

    void setThreads(int threads) {
        this->threads = threads;
    }

    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }
    
    const vector<shared_ptr<vector<uint32_t>>>& getResults() const {
        return results;
    }
    
    const vector<uint32_t>& getNNs(size_t index) const {
        return *results[index];
    }
    
    void execute();
    
    void print(ostream& out);
};

}
}
