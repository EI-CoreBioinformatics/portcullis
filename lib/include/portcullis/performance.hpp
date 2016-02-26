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
using std::stringstream;

#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

namespace portcullis {  
    
class Performance {
private:
    uint32_t tp;
    uint32_t tn;
    uint32_t fp;
    uint32_t fn;    
public:
    Performance(uint32_t tp, uint32_t tn, uint32_t fp, uint32_t fn) : tp(tp), tn(tn), fp(fp), fn(fn) {}
        
    string toString() {
        stringstream ss;
        ss << std::fixed << std::setprecision(2);
        ss << "TP:" << tp << " TN:" << tn << " FP:" << fp << " FN:" << fn 
           << " REC:" << getRecall() << " PRC:" << getPrecision() << " F1:" << getF1Score();
        return ss.str();
    }
    
    inline uint32_t getAllPositive() {
        return tp + fp;
    }

    inline uint32_t getAllNegative() {
        return tn + fn;
    }
    
    inline uint32_t getAllTrue() {
        return tn + tp;
    }
    
    inline uint32_t getAllFalse() {
        return fn + fp;
    }
        
    inline uint32_t getAll() {
        return tp + tn + fp + fn;
    }
        
    inline double getPrecision() {
        return 100.0 * (double)tp / (double)(tp + fp);
    }
    
    inline double getRecall() {
        return 100.0 * (double)tp / (double)(tp + fn);
    }
    
    // Same as recall
    inline double getSensitivity() {
        return 100.0 * (double)tp / (double)(tp + fn);
    }
    
    inline double getSpecificity() {
        return 100.0 * (double)tn / (double)(fp + tn);
    }
    
    inline double getNPV() {
        return 100.0 * (double)tn / (double)(tn + fn);
    }
    
    inline double getFallOut() {
        return 100.0 * (double)fp / (double)(fp + tn);
    }
    
    inline double getFDR() {
        return 100.0 * (double)fp / (double)(fp + tp);
    }
    
    inline double getFNR() {
        return 100.0 * (double)fn / (double)(fn + tp);
    }
    
    inline double getAccuracy() {
        return 100.0 * (double)getAllTrue() / (double)getAll();
    }
    
    inline double getF1Score() {
        return 100.0 * (double)(2.0 * tp) / (double)(2.0 * tp + fp + fn);
    }

    inline double getMCC() {
        return 100.0 * (double)(tp * tn - fp * fn) / (std::sqrt((double)(tp + fp)*(tp+fn)*(tn+fp)*(tn+fn)));
    }
    
    inline double getInformedness() {
        return getSensitivity() + getSpecificity() - 1.0;
    }
    
    inline double getMarkedness() {
        return getPrecision() + getNPV() - 1.0;
    }
    
    static void loadGenuine(path& genuineFile, vector<bool>& results) {
    
        // Load reference data    
        ifstream refs(genuineFile.string());
        string line;
        uint32_t lineNb = 0;
        while (std::getline(refs, line)) {
            std::istringstream iss(line);
            bool res;
            iss >> res;
            results.push_back(res);
        }
        refs.close();
    }
};

}