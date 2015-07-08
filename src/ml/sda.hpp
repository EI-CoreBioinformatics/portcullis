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
#include <math.h>
using namespace std;

#include "da.hpp"
#include "logistic_regression.hpp"
#include "hidden_layer.hpp"

namespace portcullis {
    namespace ml {

        class StackedDenoisingAutoencoder {
        public:
            int N;
            int n_ins;
            int *hidden_layer_sizes;
            int n_outs;
            int n_layers;
            HiddenLayer **sigmoid_layers;
            DenoisingAutoencoder **dA_layers;
            LogisticRegression *log_layer;
            StackedDenoisingAutoencoder(int, int, int*, int, int);
            ~StackedDenoisingAutoencoder();
            void pretrain(int*, double, double, int);
            void finetune(int*, int*, double, int);
            void predict(int*, double*);
        };

    }
}


