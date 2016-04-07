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

#include <boost/algorithm/string.hpp>
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
    
    static string shortHeader() {
        return "TP\tTN\tFP\tFN\tREC\tPRC\tF1";
    }
    
    static string longHeader() {
        return "TP\tTN\tFP\tFN\tSEN\tSPC\tPPV\tNPV\tF1\tACC\tINF\tMRK";
    }
    
    static string to_2dp_string(const double v) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(2) << v;
        return out.str();
    }
    
    
    string toShortString() const {
        vector<string> parts;
        parts.push_back(std::to_string(tp));
        parts.push_back(std::to_string(tn)); 
        parts.push_back(std::to_string(fp));
        parts.push_back(std::to_string(fn));
        parts.push_back(to_2dp_string(getRecall()));
        parts.push_back(to_2dp_string(getPrecision()));
        parts.push_back(to_2dp_string(getF1Score()));
        return boost::algorithm::join(parts, "\t");
    }
    
    string toLongString() const {
        vector<string> parts;
        parts.push_back(std::to_string(tp));
        parts.push_back(std::to_string(tn)); 
        parts.push_back(std::to_string(fp));
        parts.push_back(std::to_string(fn));
        parts.push_back(to_2dp_string(getSensitivity()));
        parts.push_back(to_2dp_string(getSpecificity()));
        parts.push_back(to_2dp_string(getPrecision()));
        parts.push_back(to_2dp_string(getNPV()));
        parts.push_back(to_2dp_string(getF1Score()));
        parts.push_back(to_2dp_string(getAccuracy()));
        parts.push_back(to_2dp_string(getInformedness()));
        parts.push_back(to_2dp_string(getMarkedness()));
        return boost::algorithm::join(parts, "\t");
    }
    
    
    inline uint32_t getAllPositive() const {
        return tp + fp;
    }

    inline uint32_t getAllNegative() const {
        return tn + fn;
    }
    
    inline uint32_t getAllTrue() const {
        return tn + tp;
    }
    
    inline uint32_t getAllFalse() const {
        return fn + fp;
    }
        
    inline uint32_t getAll() const {
        return tp + tn + fp + fn;
    }
        
    inline double getPrecision() const {
        return 100.0 * (double)tp / (double)(tp + fp);
    }
    
    inline double getRecall() const {
        return 100.0 * (double)tp / (double)(tp + fn);
    }
    
    // Same as recall
    inline double getSensitivity() const {
        return 100.0 * (double)tp / (double)(tp + fn);
    }
    
    inline double getSpecificity() const {
        return 100.0 * (double)tn / (double)(fp + tn);
    }
    
    inline double getNPV() const {
        return 100.0 * (double)tn / (double)(tn + fn);
    }
    
    inline double getFallOut() const {
        return 100.0 * (double)fp / (double)(fp + tn);
    }
    
    inline double getFDR() const {
        return 100.0 * (double)fp / (double)(fp + tp);
    }
    
    inline double getFNR() const {
        return 100.0 * (double)fn / (double)(fn + tp);
    }
    
    inline double getAccuracy() const {
        return 100.0 * (double)getAllTrue() / (double)getAll();
    }
    
    inline double getF1Score() const {
        return 100.0 * (double)(2.0 * tp) / (double)(2.0 * tp + fp + fn);
    }

    inline double getMCC() const {
        return 100.0 * (double)(tp * tn - fp * fn) / (std::sqrt((double)(tp + fp)*(tp+fn)*(tn+fp)*(tn+fn)));
    }
    
    inline double getInformedness() const {
        return getSensitivity() + getSpecificity() - 100.0;
    }
    
    inline double getMarkedness() const {
        return getPrecision() + getNPV() - 100.0;
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