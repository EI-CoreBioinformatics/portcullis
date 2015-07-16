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

#include <vector>
using std::vector;

#include "math.hpp"

#include "da.hpp"

portcullis::ml::DenoisingAutoencoder::DenoisingAutoencoder(int size, int n_v, int n_h, double **w, const vector<double>& hb, const vector<double>& vb) {
    N = size;
    n_visible = n_v;
    n_hidden = n_h;

    if (w == NULL) {
        W = new double*[n_hidden];
        for (int i = 0; i < n_hidden; i++) W[i] = new double[n_visible];
        double a = 1.0 / n_visible;

        for (int i = 0; i < n_hidden; i++) {
            for (int j = 0; j < n_visible; j++) {
                W[i][j] = uniform(-a, a);
            }
        }
    } else {
        W = w;
    }

    if (hb == NULL) {
        hbias = new double[n_hidden];
        for (int i = 0; i < n_hidden; i++) hbias[i] = 0;
    } else {
        hbias = hb;
    }

    if (vb == NULL) {
        vbias = new double[n_visible];
        for (int i = 0; i < n_visible; i++) vbias[i] = 0;
    } else {
        vbias = vb;
    }
}

portcullis::ml::DenoisingAutoencoder::~DenoisingAutoencoder() {
    delete[] vbias;
}

void portcullis::ml::DenoisingAutoencoder::corrupt_input(const vector<double>& in, vector<double>& out, double p) {
    for (int i = 0; i < n_visible; i++) {
        if (in[i] == 0) {
            out[i] = 0;
        } else {
            out[i] = binomial(1, p);
        }
    }
}

void portcullis::ml::DenoisingAutoencoder::encode(const vector<double>& in, vector<double>& out) {
    for (int i = 0; i < n_hidden; i++) {
        double val = 0;
        for (int j = 0; j < n_visible; j++) {
            val += W[i][j] * in[j];
        }
        val += hbias[i];
        out[i] = sigmoid(val);
    }
}

void portcullis::ml::DenoisingAutoencoder::decode(const vector<double>& in, vector<double>& out) {
    for (int i = 0; i < n_visible; i++) {
        double val = 0;
        for (int j = 0; j < n_hidden; j++) {
            val += W[j][i] * in[j];
        }
        val += vbias[i];
        out[i] = sigmoid(val);
    }
}

void portcullis::ml::DenoisingAutoencoder::train(const vector<double>& x, const double learning_rate, const double corruption_level) {
    vector<double> tilde_x(n_visible);
    double *y = new double[n_hidden];
    double *z = new double[n_visible];

    double *L_vbias = new double[n_visible];
    double *L_hbias = new double[n_hidden];

    const double p = 1 - corruption_level;

    corrupt_input(x, tilde_x, p);
    encode(tilde_x, y);
    decode(y, z);

    // vbias
    for (int i = 0; i < n_visible; i++) {
        L_vbias[i] = x[i] - z[i];
        vbias[i] += learning_rate * L_vbias[i] / N;
    }

    // hbias
    for (int i = 0; i < n_hidden; i++) {
        L_hbias[i] = 0;
        for (int j = 0; j < n_visible; j++) {
            L_hbias[i] += W[i][j] * L_vbias[j];
        }
        L_hbias[i] *= y[i] * (1 - y[i]);

        hbias[i] += learning_rate * L_hbias[i] / N;
    }

    // W
    for (int i = 0; i < n_hidden; i++) {
        for (int j = 0; j < n_visible; j++) {
            W[i][j] += learning_rate * (L_hbias[i] * tilde_x[j] + L_vbias[j] * y[i]) / N;
        }
    }

    delete[] L_hbias;
    delete[] L_vbias;
    delete[] z;
    delete[] y;
}

void portcullis::ml::DenoisingAutoencoder::reconstruct(const vector<double>& x, vector<double>& z) {    
    vector<double> y(n_hidden);
    encode(x, y);
    decode(y, z);
}
