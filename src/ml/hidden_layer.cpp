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
#include <math.h>
using namespace std;

#include "math.hpp"

#include "hidden_layer.hpp"

portcullis::ml::HiddenLayer::HiddenLayer(int size, int in, int out, double **w, double *bp) {
    N = size;
    n_in = in;
    n_out = out;

    if (w == NULL) {
        W = new double*[n_out];
        for (int i = 0; i < n_out; i++) W[i] = new double[n_in];
        double a = 1.0 / n_in;

        for (int i = 0; i < n_out; i++) {
            for (int j = 0; j < n_in; j++) {
                W[i][j] = uniform(-a, a);
            }
        }
    } else {
        W = w;
    }

    if (bp == NULL) {
        b = new double[n_out];
    } else {
        b = bp;
    }
}

portcullis::ml::HiddenLayer::~HiddenLayer() {
    for (int i = 0; i < n_out; i++) delete W[i];
    delete[] W;
    delete[] b;
}

double portcullis::ml::HiddenLayer::output(int *input, const vector<double>& w, double b) {
    double linear_output = 0.0;
    for (int j = 0; j < n_in; j++) {
        linear_output += w[j] * input[j];
    }
    linear_output += b;
    return sigmoid(linear_output);
}

void portcullis::ml::HiddenLayer::sample_h_given_v(int *input, int *sample) {
    for (int i = 0; i < n_out; i++) {
        sample[i] = binomial(1, output(input, W[i], b[i]));
    }
}
