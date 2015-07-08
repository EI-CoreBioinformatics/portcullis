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

#include <stdlib.h>
#include <math.h>

double uniform(double min, double max) {
    return rand() / (RAND_MAX + 1.0) * (max - min) + min;
}

int binomial(int n, double p) {
    if (p < 0.0 || p > 1.0) return 0;

    int c = 0;
    double r = 0.0;

    for (int i = 0; i < n; i++) {
        r = rand() / (RAND_MAX + 1.0);
        if (r < p) c++;
    }

    return c;
}

double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}



/*
void test_sda() {
  srand(0);

  double pretrain_lr = 0.1;
  double corruption_level = 0.3;
  int pretraining_epochs = 1000;
  double finetune_lr = 0.1;
  int finetune_epochs = 500;

  int train_N = 10;
  int test_N = 4;
  int n_ins = 28;
  int n_outs = 2;
  int hidden_layer_sizes[] = {15, 15};
  int n_layers = sizeof(hidden_layer_sizes) / sizeof(hidden_layer_sizes[0]);

  // training data
  int train_X[10][28] = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1}
  };

  int train_Y[10][2] = {
    {1, 0},
    {1, 0},
    {1, 0},
    {1, 0},
    {1, 0},
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1},
    {0, 1}
  };

  // construct SdA
  SdA sda(train_N, n_ins, hidden_layer_sizes, n_outs, n_layers);

  // pretrain
  sda.pretrain(*train_X, pretrain_lr, corruption_level, pretraining_epochs);

  // finetune
  sda.finetune(*train_X, *train_Y, finetune_lr, finetune_epochs);


  // test data
  int test_X[4][28] = {
    {1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}
  };

  double test_Y[4][28];

  // test
  for(int i=0; i<test_N; i++) {
    sda.predict(test_X[i], test_Y[i]);
    for(int j=0; j<n_outs; j++) {
      printf("%.5f ", test_Y[i][j]);
    }
    cout << endl;
  }
  
}


int main() {
  test_sda();
  return 0;
}*/
