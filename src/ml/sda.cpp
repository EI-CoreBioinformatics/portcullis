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


#include "sda.hpp"

portcullis::ml::StackedDenoisingAutoencoder::StackedDenoisingAutoencoder(int size, int n_i, int *hls, int n_o, int n_l) {
    int input_size;

    N = size;
    n_ins = n_i;
    hidden_layer_sizes = hls;
    n_outs = n_o;
    n_layers = n_l;

    sigmoid_layers = new HiddenLayer*[n_layers];
    dA_layers = new DenoisingAutoencoder*[n_layers];

    // construct multi-layer
    for (int i = 0; i < n_layers; i++) {
        if (i == 0) {
            input_size = n_ins;
        } else {
            input_size = hidden_layer_sizes[i - 1];
        }

        // construct sigmoid_layer
        sigmoid_layers[i] = new HiddenLayer(N, input_size, hidden_layer_sizes[i], NULL, NULL);

        // construct dA_layer
        dA_layers[i] = new DenoisingAutoencoder(N, input_size, hidden_layer_sizes[i],\
                          sigmoid_layers[i]->W, sigmoid_layers[i]->b, NULL);
    }

    // layer for output using LogisticRegression
    log_layer = new LogisticRegression(N, hidden_layer_sizes[n_layers - 1], n_outs);
}

portcullis::ml::StackedDenoisingAutoencoder::~StackedDenoisingAutoencoder() {
    delete log_layer;

    for (int i = 0; i < n_layers; i++) {
        delete sigmoid_layers[i];
        delete dA_layers[i];
    }
    delete[] sigmoid_layers;
    delete[] dA_layers;
}

void portcullis::ml::StackedDenoisingAutoencoder::pretrain(int *input, const double learning_rate, double corruption_level, int epochs) {
    int *layer_input;
    int prev_layer_input_size;
    int *prev_layer_input;

    int *train_X = new int[n_ins];

    for (int i = 0; i < n_layers; i++) { // layer-wise

        for (int epoch = 0; epoch < epochs; epoch++) { // training epochs

            for (int n = 0; n < N; n++) { // input x1...xN
                // initial input
                for (int m = 0; m < n_ins; m++) train_X[m] = input[n * n_ins + m];

                // layer input
                for (int l = 0; l <= i; l++) {

                    if (l == 0) {
                        layer_input = new int[n_ins];
                        for (int j = 0; j < n_ins; j++) layer_input[j] = train_X[j];
                    } else {
                        if (l == 1) prev_layer_input_size = n_ins;
                        else prev_layer_input_size = hidden_layer_sizes[l - 2];

                        prev_layer_input = new int[prev_layer_input_size];
                        for (int j = 0; j < prev_layer_input_size; j++) prev_layer_input[j] = layer_input[j];
                        delete[] layer_input;

                        layer_input = new int[hidden_layer_sizes[l - 1]];

                        sigmoid_layers[l - 1]->sample_h_given_v(prev_layer_input, layer_input);
                        delete[] prev_layer_input;
                    }
                }

                dA_layers[i]->train(layer_input, learning_rate, corruption_level);

            }
        }
    }

    delete[] train_X;
    delete[] layer_input;
}

void portcullis::ml::StackedDenoisingAutoencoder::finetune(int *input, int *label, double lr, int epochs) {
    int *layer_input;
    int prev_layer_input_size;
    int *prev_layer_input;

    int *train_X = new int[n_ins];
    int *train_Y = new int[n_outs];

    for (int epoch = 0; epoch < epochs; epoch++) {
        for (int n = 0; n < N; n++) { // input x1...xN
            // initial input
            for (int m = 0; m < n_ins; m++) train_X[m] = input[n * n_ins + m];
            for (int m = 0; m < n_outs; m++) train_Y[m] = label[n * n_outs + m];

            // layer input
            for (int i = 0; i < n_layers; i++) {
                if (i == 0) {
                    prev_layer_input = new int[n_ins];
                    for (int j = 0; j < n_ins; j++) prev_layer_input[j] = train_X[j];
                } else {
                    prev_layer_input = new int[hidden_layer_sizes[i - 1]];
                    for (int j = 0; j < hidden_layer_sizes[i - 1]; j++) prev_layer_input[j] = layer_input[j];
                    delete[] layer_input;
                }


                layer_input = new int[hidden_layer_sizes[i]];
                sigmoid_layers[i]->sample_h_given_v(prev_layer_input, layer_input);
                delete[] prev_layer_input;
            }

            log_layer->train(layer_input, train_Y, lr);
        }
        // lr *= 0.95;
    }

    delete[] layer_input;
    delete[] train_X;
    delete[] train_Y;
}

void portcullis::ml::StackedDenoisingAutoencoder::predict(int *x, double *y) {
    double *layer_input;
    int prev_layer_input_size;
    double *prev_layer_input;

    double linear_output;

    prev_layer_input = new double[n_ins];
    for (int j = 0; j < n_ins; j++) prev_layer_input[j] = x[j];

    // layer activation
    for (int i = 0; i < n_layers; i++) {
        layer_input = new double[sigmoid_layers[i]->n_out];

        for (int k = 0; k < sigmoid_layers[i]->n_out; k++) {
            linear_output = 0.0;

            for (int j = 0; j < sigmoid_layers[i]->n_in; j++) {
                linear_output += sigmoid_layers[i]->W[k][j] * prev_layer_input[j];
            }
            linear_output += sigmoid_layers[i]->b[k];
            layer_input[k] = sigmoid(linear_output);
        }
        delete[] prev_layer_input;

        if (i < n_layers - 1) {
            prev_layer_input = new double[sigmoid_layers[i]->n_out];
            for (int j = 0; j < sigmoid_layers[i]->n_out; j++) prev_layer_input[j] = layer_input[j];
            delete[] layer_input;
        }
    }

    for (int i = 0; i < log_layer->n_out; i++) {
        y[i] = 0;
        for (int j = 0; j < log_layer->n_in; j++) {
            y[i] += log_layer->W[i][j] * layer_input[j];
        }
        y[i] += log_layer->b[i];
    }

    log_layer->softmax(y);


    delete[] layer_input;
}
