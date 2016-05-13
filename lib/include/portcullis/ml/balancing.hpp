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
    
typedef boost::error_info<struct BalancingError,string> BalancingErrorInfo;
struct BalacingException: virtual boost::exception, virtual std::exception { };

class UnderSampler : UnbalancedDataset {
    
};

/**
 * Object to perform classification on balanced ensembled selected from random sampling.
 * 
 * It is based on the idea presented in the paper "Exploratory Undersampling 
 * Class-Imbalance Learning" by Liu et al.    
 */
class EasyEnsemble : public UnderSampler {
private:
     uint16_t nb_subsets;
     
     
        
public:

    EasyEnsemble() {
        
    }
};
}
}