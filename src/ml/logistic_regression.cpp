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

#include "logistic_regression.hpp"

portcullis::ml::LogisticRegression::LogisticRegression(int size, int in, int out) {
    N = size;
    n_in = in;
    n_out = out;

    W = new double*[n_out];
    for (int i = 0; i < n_out; i++) W[i] = new double[n_in];
    b = new double[n_out];

    for (int i = 0; i < n_out; i++) {
        for (int j = 0; j < n_in; j++) {
            W[i][j] = 0;
        }
        b[i] = 0;
    }
}

portcullis::ml::LogisticRegression::~LogisticRegression() {
    for (int i = 0; i < n_out; i++) delete[] W[i];
    delete[] W;
    delete[] b;
}

void portcullis::ml::LogisticRegression::train(int *x, int *y, const double learning_rate) {
    double *p_y_given_x = new double[n_out];
    double *dy = new double[n_out];

    for (int i = 0; i < n_out; i++) {
        p_y_given_x[i] = 0;
        for (int j = 0; j < n_in; j++) {
            p_y_given_x[i] += W[i][j] * x[j];
        }
        p_y_given_x[i] += b[i];
    }
    softmax(p_y_given_x);

    for (int i = 0; i < n_out; i++) {
        dy[i] = y[i] - p_y_given_x[i];

        for (int j = 0; j < n_in; j++) {
            W[i][j] += learning_rate * dy[i] * x[j] / N;
        }

        b[i] += learning_rate * dy[i] / N;
    }

    delete[] p_y_given_x;
    delete[] dy;
}

void portcullis::ml::LogisticRegression::softmax(vector<double>& x) {
    double max = 0.0;
    double sum = 0.0;

    for (int i = 0; i < n_out; i++) if (max < x[i]) max = x[i];
    for (int i = 0; i < n_out; i++) {
        x[i] = exp(x[i] - max);
        sum += x[i];
    }

    for (int i = 0; i < n_out; i++) x[i] /= sum;
}

void portcullis::ml::LogisticRegression::predict(int *x, vector<double>& y) {
    for (int i = 0; i < n_out; i++) {
        y[i] = 0;
        for (int j = 0; j < n_in; j++) {
            y[i] += W[i][j] * x[j];
        }
        y[i] += b[i];
    }

    softmax(y);
}
