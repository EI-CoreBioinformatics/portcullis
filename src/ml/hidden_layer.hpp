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

        class HiddenLayer {
        public:
            int N;
            int n_in;
            int n_out;
            double **W;
            double *b;
            HiddenLayer(int, int, int, double**, double*);
            ~HiddenLayer();
            double output(int*, const vector<double>&, double);
            void sample_h_given_v(int*, int*);
        };

    }
}