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
#include <thread>
#include <vector>
using std::string;
using std::thread;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

#include <ranger/Data.h>

namespace portcullis {
namespace ml {

typedef boost::error_info<struct ENNError,string> ENNErrorInfo;
struct ENNException: virtual boost::exception, virtual std::exception { };
    
const uint16_t ENN_THRESHOLD = 5;

/**
 * An parallel implementation of Wilson's Edited nearest neighbour algorithm.
 * Uses KNN to identify nearest neighbours for each sample, then marks entries
 * to be removed or to be kept if more than half of the KNNs come from a different
 * class.
 * KNN logic originally derived from OpenCV, modified to work with ranger data
 * types 
 */
class ENN {
protected:
    uint16_t k;
    uint16_t threshold;
    uint16_t threads;
    bool verbose;
    
    double* data;
    size_t rows;
    size_t cols;
    
    vector<bool> labels;    
    
public:    
    
    ENN(uint16_t defaultK, uint16_t _threads, double* _data, size_t _rows, size_t _cols, vector<bool>& _labels);
    
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
    
    void setThreshold(uint16_t threshold) {
        this->threshold = threshold;
    }
    
    uint16_t getThreshold() const {
        return threshold;
    }

    
    uint32_t execute( vector<bool>& results ) const;
};

}
}
