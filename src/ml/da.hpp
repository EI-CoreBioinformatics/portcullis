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

namespace portcullis {
    namespace ml {

        class DenoisingAutoencoder {
        public:
            int N;
            int n_visible;
            int n_hidden;
            double **W;
            double *hbias;
            double *vbias;
            DenoisingAutoencoder(int, int, int, double**, double*, double*);
            ~DenoisingAutoencoder();
            void corrupt_input(const vector<double>&, int*, double);
            void encode(const vector<double>&, vector<double>&);
            void decode(const vector<double>&, vector<double>&);
            void train(int*, double, double);
            void reconstruct(const vector<double>&, vector<double>&);
        };

    }
}